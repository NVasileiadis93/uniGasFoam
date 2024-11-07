/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "uniGasMassFlowRateInletPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{
defineTypeNameAndDebug(uniGasMassFlowRateInletPatch, 0);

addToRunTimeSelectionTable
(
    uniGasGeneralBoundary,
    uniGasMassFlowRateInletPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasMassFlowRateInletPatch::uniGasMassFlowRateInletPatch
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    theta_(propsDict_.getOrDefault<scalar>("theta",1.0)),
    inletTemperature_(propsDict_.get<scalar>("inletTemperature")),
    massFlowRate_(propsDict_.get<scalar>("massFlowRate")),
    initialVelocity_(propsDict_.getOrDefault<vector>("initialVelocity",vector::zero)),
    moleFractions_(),
    patchSurfaceArea_(),
    moleFlowRate_(),
    parcelsIn_(),
    parcelsToInsert_(),
    inletNumberDensity_(),
    inletVelocity_(faces_.size(), initialVelocity_),
    previousInletVelocity_(faces_.size(), initialVelocity_)

{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    // Compute total patch area
    patchSurfaceArea_ = 0.0;
    forAll(faces_,f)
    {
        patchSurfaceArea_ += mag(mesh.faceAreas()[faces_[f]]);
    }

    if (Pstream::parRun())
    {
        reduce(patchSurfaceArea_, sumOp<scalar>());
    }

    Info << "TOTAL PATCH AREA: " << patchSurfaceArea_ << endl;

    // Get type IDs
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // Set the accumulator
    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(), 0.0);
    }

    // read in the mole fraction per specie
    const dictionary& moleFractionsDict = propsDict_.subDict("moleFractions");

    moleFractions_.clear();

    moleFractions_.setSize(typeIds_.size(), 0.0);

    forAll(moleFractions_, i)
    {
        const word& moleculeName = cloud_.typeIdList()[typeIds_[i]];
        moleFractions_[i] = moleFractionsDict.get<scalar>(moleculeName);
    }

    inletNumberDensity_.setSize(typeIds_.size());
    
    forAll(inletNumberDensity_, m)
    {
        inletNumberDensity_[m].setSize(faces_.size(), 0.0);
    }

    moleFlowRate_.setSize(typeIds_.size()); 

    parcelsIn_.setSize(typeIds_.size());
    
    parcelsToInsert_.setSize(typeIds_.size());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasMassFlowRateInletPatch::initialConfiguration()
{}


void Foam::uniGasMassFlowRateInletPatch::calculateProperties()
{}


void Foam::uniGasMassFlowRateInletPatch::controlParcelsBeforeMove()
{

    previousInletVelocity_ = inletVelocity_;
    inletVelocity_ = vector::zero;

    insertParcels
    (
        inletTemperature_,
        inletVelocity_
    );

}


void Foam::uniGasMassFlowRateInletPatch::controlParcelsBeforeCollisions()
{
}


void Foam::uniGasMassFlowRateInletPatch::controlParcelsAfterCollisions()
{

    const scalar deltaT = mesh_.time().deltaTValue();
    const scalar nParticle = cloud_.nParticle();

    vectorField newInletVelocity_(faces_.size(), Zero);
    vectorField momentum(faces_.size(), Zero);
    scalarField mass(faces_.size(), scalar(0));

    const List<DynamicList<uniGasParcel*> >& cellOccupancy = cloud_.cellOccupancy();
    const List<scalarField>& parcelIdFlux = cloud_.tracker().parcelIdFlux();

    scalar totalMass = 0;
    forAll(moleFractions_, iD)
    {
        const label typeId = typeIds_[iD];
        const scalar pMass = cloud_.constProps(typeId).mass();
        totalMass += pMass*moleFractions_[iD];
    }

    forAll(moleFractions_, iD)
    {
        inletNumberDensity_[iD] = 0.0;
        parcelsIn_[iD] = 0.0;
        moleFlowRate_[iD] = moleFractions_[iD]*(massFlowRate_ / totalMass);
    }

    forAll(cells_, c)
    {

        const label faceI = faces_[c];
        const label cellI = cells_[c];

        forAll(moleFractions_, iD)
        {
            scalar CWF = cloud_.cellWF(cellI);
            scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);
            parcelsIn_[iD] += moleFlowRate_[iD]*deltaT*(mag(mesh_.faceAreas()[faceI])/patchSurfaceArea_)/(nParticle*CWF*RWF) + parcelIdFlux[iD][faceI]/(CWF*RWF);
        }

        const List<uniGasParcel*>& parcelsInCell = cellOccupancy[cellI];

        // compute cell instantaneous numberDensity and velocity
        for (uniGasParcel* p : parcelsInCell)
        {
            const label iD = p->typeId();
            
            const scalar pMass = nParticle*cloud_.constProps(iD).mass();

            scalar CWF = cloud_.cellWF(cellI);
            scalar RWF = cloud_.axiRWF(p->position());

            inletNumberDensity_[iD][c] += 1;
            momentum[c] += pMass*CWF*RWF*p->U();
            mass[c] += pMass*CWF*RWF;

        }

        if (mass[c] > VSMALL)
        {
            newInletVelocity_[c] = momentum[c]/mass[c];
        }

        inletVelocity_[c] = theta_*newInletVelocity_[c] + (1.0-theta_)*previousInletVelocity_[c];

        const vector& sF = mesh_.faceAreas()[faceI];
        const scalar fA = mag(sF);    

        if ((inletVelocity_[c] & -sF/fA) < 0)
        {
            inletVelocity_[c] = previousInletVelocity_[c];
        }

        forAll(moleFractions_, iD)
        {

            scalar CWF = cloud_.cellWF(cellI);
            scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);
            inletNumberDensity_[iD][c] = inletNumberDensity_[iD][c]*nParticle*CWF*RWF/mesh_.V()[cellI] ;

        }

    }

    parcelsToInsert_ = 0;
    forAll(cells_, c)
    {

        const label faceI = faces_[c];
        const label cellI = cells_[c];

        forAll(moleFractions_, iD)
        {

            const label typeId = typeIds_[iD];

            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            scalar pMass = cloud_.constProps(typeId).mass();

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    inletTemperature_,
                    pMass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (inletVelocity_[c] & -sF/fA )/mostProbableSpeed;

            scalar CWF = cloud_.cellWF(cellI);
            scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);

            parcelsToInsert_ += 
            (
                fA*inletNumberDensity_[iD][c]*deltaT*mostProbableSpeed
               *(
                    exp(-sqr(sCosTheta))
                  + sqrt(pi)*sCosTheta*(1 + erf(sCosTheta))
                )
            )
           /(2.0*sqrt(pi)*nParticle*CWF*RWF);
           
        }

    }

    if (Pstream::parRun())
    {
        reduce(parcelsIn_, sumOp<scalarField>());
        reduce(parcelsToInsert_, sumOp<scalarField>());
    }

    forAll(cells_, c)
    {

        forAll(moleFractions_, iD)
        {
            inletNumberDensity_[iD][c] = inletNumberDensity_[iD][c]*(parcelsIn_[iD]/parcelsToInsert_[iD]);
        }
    
    }

    computeParcelsToInsert
    (
        inletNumberDensity_,
        inletTemperature_,
        inletVelocity_
    );

}


void Foam::uniGasMassFlowRateInletPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uniGasMassFlowRateInletPatch::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uniGasGeneralBoundary::updateProperties(dict);
}

// ************************************************************************* //