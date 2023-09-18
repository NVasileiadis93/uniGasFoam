/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "uspWangPressureInlet.H"
#include "uspCloud.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{
defineTypeNameAndDebug(uspWangPressureInlet, 0);

addToRunTimeSelectionTable
(
    uspGeneralBoundary,
    uspWangPressureInlet,
    dictionary
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspWangPressureInlet::uspWangPressureInlet
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    moleFractions_(),
    inletPressure_(),
    inletTemperature_(),
    n_(),
    cellVolume_(faces_.size(), scalar(00)),
    inletVelocity_(faces_.size(), Zero),
    totalMomentum_(faces_.size(), Zero),
    totalMass_(faces_.size(), scalar(0)),
    nTotalParcels_(faces_.size(), scalar(0)),
    nTimeSteps_(scalar(0)),
    mcc_(faces_.size(), scalar(0.0)),
    UMean_(faces_.size(), Zero),
    UCollected_(faces_.size(), Zero),
    velSqrMean_(faces_.size(), Zero),
    velMeanSqr_(faces_.size(), Zero),
    totalVelSqrMean_(faces_.size(), Zero),
    totalVelMeanSqr_(faces_.size(), Zero)
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();

    // Calculate required number density at inlet boundary
    // equation 32, Liou and Fang (2000)
    n_ = inletPressure_ / (physicoChemical::k.value()*inletTemperature_);

     // get volume of each boundary cell
    forAll(cellVolume_, c)
    {
        cellVolume_[c] = mesh_.cellVolumes()[cells_[c]];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspWangPressureInlet::initialConfiguration()
{}


void Foam::uspWangPressureInlet::calculateProperties()
{}


void Foam::uspWangPressureInlet::controlParcelsBeforeMove()
{
    insertParcels
    (
        inletTemperature_,
        inletVelocity_
    );
}


void Foam::uspWangPressureInlet::controlParcelsBeforeCollisions()
{}


void Foam::uspWangPressureInlet::controlParcelsAfterCollisions()
{
    nTimeSteps_ += 1.0;

    scalar molecularMass = 0.0;
    scalar molarconstantPressureSpecificHeat = 0.0;
    scalar molarconstantVolumeSpecificHeat = 0.0;

    forAll(moleFractions_, iD)
    {
        const label typeId_ = typeIds_[iD];

        const auto& constProps = cloud_.constProps(typeId_);

        molecularMass += constProps.mass()*moleFractions_[iD];
        molarconstantPressureSpecificHeat +=
            (5.0 + constProps.rotationalDoF())*moleFractions_[iD];
        molarconstantVolumeSpecificHeat +=
            (3.0 + constProps.rotationalDoF())*moleFractions_[iD];
    }

     // R = k/m
    const scalar gasConstant = physicoChemical::k.value()/molecularMass;

    const scalar gamma =
        molarconstantPressureSpecificHeat
       /molarconstantVolumeSpecificHeat;

    vectorField momentum(faces_.size(),Zero);
    vectorField UCollected(faces_.size(), Zero);
    scalarField mass(faces_.size(), Zero);
    scalarField nParcels(faces_.size(), Zero);
    scalarField mcc(faces_.size(), Zero);
    scalarField translationalTemperature(faces_.size(), Zero);
    scalarField numberDensity(faces_.size(), Zero);
    scalarField massDensity(faces_.size(), Zero);
    scalarField pressure(faces_.size(), Zero);
    scalarField speedOfSound(faces_.size(), Zero);
    scalarField velocityCorrection(faces_.size(), Zero);
    vectorField velSqrMean(faces_.size(), Zero);
    vectorField velMeanSqr(faces_.size(), Zero);

    const auto& cellOccupancy = cloud_.cellOccupancy();

    const scalar nParticle = cloud_.nParticle();

    forAll(cells_, c)
    {

        const label cellI = cells_[c];

        for (uspParcel* p : cellOccupancy[cells_[c]])
        {
            label iD = typeIds_.find(p->typeId());

            if (iD != -1)
            {
                const auto& constProps = cloud_.constProps(p->typeId());

                scalar CWF = cloud_.cellWF(cellI);
                scalar RWF = cloud_.axiRWF(p->position());

                momentum[c] += nParticle*CWF*RWF*constProps.mass()*p->U();
                mass[c] += nParticle*CWF*RWF*constProps.mass();

                nParcels[c] += 1.0;
                UCollected[c] += p->U();
                velSqrMean[c] += cmptSqr(p->U());
                velMeanSqr[c] += p->U();
            }
        }

        nTotalParcels_[c] += nParcels[c];

        totalMomentum_[c] += momentum[c];

        totalMass_[c] += mass[c];

        mcc_[c] += mcc[c];

        UCollected_[c] += UCollected[c];

        massDensity[c] = totalMass_[c]/(cellVolume_[c]*nTimeSteps_);

        numberDensity[c] = massDensity[c]/molecularMass;

        if (nTotalParcels_[c] > 1)
        {
            velSqrMean_[c] += velSqrMean[c];
            velMeanSqr_[c] += velMeanSqr[c];

            totalVelSqrMean_[c] = velSqrMean_[c]/nTotalParcels_[c];
            totalVelMeanSqr_[c] = cmptSqr(velMeanSqr_[c]/nTotalParcels_[c]);

            UMean_[c] = UCollected_[c]/nTotalParcels_[c];

            translationalTemperature[c] =
                (0.5*molecularMass)*(2.0/(3.0*physicoChemical::k.value()))
               *(
                    cmptSum(totalVelSqrMean_[c])
                   -cmptSum(totalVelMeanSqr_[c])
                );

            if (translationalTemperature[c] < VSMALL)
            {
                translationalTemperature[c] = 300.00;
            }

            pressure[c] =
                numberDensity[c]*physicoChemical::k.value()
               *translationalTemperature[c];

            speedOfSound[c] =
                sqrt(gamma*gasConstant*translationalTemperature[c]);

            label faceI = faces_[c];
            vector n = mesh_.faceAreas()[faceI];
            n.normalise();

            inletVelocity_[c] = totalMomentum_[c]/totalMass_[c];

            if (nTimeSteps_ > 100)
            {
                // Velocity correction for each boundary cellI

                velocityCorrection[c] =
                    (pressure[c] - inletPressure_)
                   /(massDensity[c]*speedOfSound[c]);
                inletVelocity_[c] += velocityCorrection[c]*n;
            }
        }
    }

    // Compute number of parcels to insert
    computeParcelsToInsert
    (
        inletTemperature_,
        inletVelocity_,
        n_,
        moleFractions_
    );
}


void Foam::uspWangPressureInlet::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uspWangPressureInlet::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uspGeneralBoundary::updateProperties(dict);

    setProperties();
}


void Foam::uspWangPressureInlet::setProperties()
{
    inletPressure_ = propsDict_.get<scalar>("inletPressure");

    inletTemperature_ = propsDict_.get<scalar>("inletTemperature");

    // Read in the type ids
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // Read in the mole fraction per specie

    const dictionary& moleFractionsDict(propsDict_.subDict("moleFractions"));

    moleFractions_.clear();

    moleFractions_.setSize(typeIds_.size(), Zero);

    forAll(moleFractions_, i)
    {
        const word& moleculeName = cloud_.typeIdList()[typeIds_[i]];
        moleFractions_[i] = moleFractionsDict.get<scalar>(moleculeName);
    }

    // Set the accumulator

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(), 0.0);
    }
}



// ************************************************************************* //
