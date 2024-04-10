/*---------------------------------------------------------------------------* \
  =========                 |
  \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration     |
    \  /    A nd           | www.openfoam.com
     \/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "uspMassFlowRateInletPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{
defineTypeNameAndDebug(uspMassFlowRateInletPatch, 0);

addToRunTimeSelectionTable
(
    uspGeneralBoundary,
    uspMassFlowRateInletPatch,
    dictionary
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspMassFlowRateInletPatch::uspMassFlowRateInletPatch
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    moleFractions_(),
    massFlowRate_(propsDict_.get<scalar>("massFlowRate")),
    inletTemperature_(propsDict_.get<scalar>("inletTemperature")),
    initialVelocity_(propsDict_.get<vector>("initialVelocity")),
    theta_(propsDict_.getOrDefault<scalar>("theta",1.0)),
    moleFlowRate_(),
    parcelsIn_(),
    parcelsOut_(),
    parcelsToInsert_(),
    inletNumberDensity_(),
    inletVelocity_(faces_.size(), Zero),
    previousInletVelocity_(faces_.size(), Zero)
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

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
    
    parcelsOut_.setSize(typeIds_.size());
    
    parcelsToInsert_.setSize(typeIds_.size());

    //sorts issues with the velocity pointing out of the mesh
    inletVelocity_ = initialVelocity_;

    previousInletVelocity_ = initialVelocity_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspMassFlowRateInletPatch::initialConfiguration()
{}


void Foam::uspMassFlowRateInletPatch::calculateProperties()
{}


void Foam::uspMassFlowRateInletPatch::controlParcelsBeforeMove()
{

    previousInletVelocity_ = inletVelocity_;
    inletVelocity_ = vector::zero;

    insertParcels
    (
        inletTemperature_,
        inletVelocity_
    );

}


void Foam::uspMassFlowRateInletPatch::controlParcelsBeforeCollisions()
{}


void Foam::uspMassFlowRateInletPatch::controlParcelsAfterCollisions()
{

    const scalar deltaT = mesh_.time().deltaTValue();
    const scalar nParticle = cloud_.nParticle();

    vectorField newInletVelocity_(faces_.size(), Zero);
    vectorField momentum(faces_.size(), Zero);
    scalarField mass(faces_.size(), scalar(0));

    const List<DynamicList<uspParcel*> >& cellOccupancy = cloud_.cellOccupancy();
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
        parcelsOut_[iD] = 0.0;
        moleFlowRate_[iD] = moleFractions_[iD]*(massFlowRate_ / totalMass);
    }

    forAll(cells_, c)
    {

        const label faceI = faces_[c];
        const label cellI = cells_[c];

        forAll(moleFractions_, iD)
        {
            parcelsOut_[iD] += parcelIdFlux[iD][faceI];
        }

        const List<uspParcel*>& parcelsInCell = cellOccupancy[cellI];

        // compute cell instantaneous numberDensity and velocity
        for (uspParcel* p : parcelsInCell)
        {
            const label iD = p->typeId();
            
            const scalar pMass = nParticle*cloud_.constProps(iD).mass();

            inletNumberDensity_[iD][c] += 1;
            momentum[c] += pMass*p->U();
            mass[c] += pMass;

        }

        if (mass[c] > VSMALL)
        {
            newInletVelocity_[c] = momentum[c]/mass[c];
            inletVelocity_[c] = theta_*newInletVelocity_[c] + (1.0-theta_)*previousInletVelocity_[c];
        }
        else
        {
            inletVelocity_[c] = previousInletVelocity_[c];
        }
        
        const vector& sF = mesh_.faceAreas()[faceI];
        const scalar fA = mag(sF);    

        if ((inletVelocity_[c] & -sF/fA) < 0)
        {
            inletVelocity_[c] = previousInletVelocity_[c];
        }

        scalar CWF = cloud_.cellWF(cellI);
        scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);
        forAll(moleFractions_, iD)
        {
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
        reduce(parcelsOut_, sumOp<scalarField>());
        reduce(parcelsToInsert_, sumOp<scalarField>());
    }

    parcelsIn_ = moleFlowRate_*deltaT/nParticle + parcelsOut_;

    
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


void Foam::uspMassFlowRateInletPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uspMassFlowRateInletPatch::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uspGeneralBoundary::updateProperties(dict);
}

// ************************************************************************* //
