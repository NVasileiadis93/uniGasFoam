/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*--------------------------------------------------------------------------*/

#include "uspLiouFangPressureOutlet.H"
#include "uspCloud.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{
defineTypeNameAndDebug(uspLiouFangPressureOutlet, 0);
addToRunTimeSelectionTable
(
    uspGeneralBoundary,
    uspLiouFangPressureOutlet,
    dictionary
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspLiouFangPressureOutlet::uspLiouFangPressureOutlet
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    outletPressure_(),
    nTimeSteps_(scalar(0.0)),
    UMean_(faces_.size(), Zero),
    outletVelocity_(faces_.size(), Zero),
    UCollected_(faces_.size(), Zero),
    totalVibrationalEnergy_(),
    vibT_(),
    vDof_(),
    nTotalParcelsSpecies_(),
    moleFractions_(),
    cellVolume_(faces_.size(), Zero),
    outletNumberDensity_(faces_.size(), Zero),
    outletMassDensity_(faces_.size(), Zero),
    outletTemperature_(faces_.size(), scalar(300.0)),
    totalMomentum_(faces_.size(), Zero),
    totalMass_(faces_.size(), Zero),
    totalMassDensity_(faces_.size(), Zero),
    totalNumberDensity_(faces_.size(), Zero),
    totalRotationalEnergy_(faces_.size(), Zero),
    totalRotationalDof_(faces_.size(), Zero),
    nTotalParcels_(faces_.size(), Zero),
    nTotalParcelsInt_(faces_.size(), Zero),
    velSqrMean_(faces_.size(), Zero),
    velMeanSqr_(faces_.size(), Zero),
    totalVelSqrMean_(faces_.size(), Zero),
    totalVelMeanSqr_(faces_.size(), Zero)
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();

    forAll(cellVolume_, c)
    {
        //get volume of each boundary cell
        cellVolume_[c] = mesh_.cellVolumes()[cells_[c]];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspLiouFangPressureOutlet::initialConfiguration()
{}


void Foam::uspLiouFangPressureOutlet::calculateProperties()
{}


void Foam::uspLiouFangPressureOutlet::controlParcelsBeforeMove()
{
    insertParcels
    (
        outletTemperature_,
        outletVelocity_
    );
}


void Foam::uspLiouFangPressureOutlet::controlParcelsBeforeCollisions()
{}


void Foam::uspLiouFangPressureOutlet::controlParcelsAfterCollisions()
{
    nTimeSteps_ += 1.0;

    scalar molecularMass = 0.0;
    scalar molarcontantPressureSpecificHeat = 0.0;
    scalar molarcontantVolumeSpecificHeat = 0.0;

    forAll(moleFractions_, iD)
    {
        const label typeId = typeIds_[iD];

        molecularMass += cloud_.constProps(typeId).mass()*moleFractions_[iD];
        molarcontantPressureSpecificHeat +=
            (5.0 + cloud_.constProps(typeId).rotationalDoF())
           *moleFractions_[iD];
        molarcontantVolumeSpecificHeat +=
            (3.0 + cloud_.constProps(typeId).rotationalDoF())
           *moleFractions_[iD];
    }

    // R = k/m
    const scalar gasConstant = physicoChemical::k.value()/molecularMass;

    const scalar gamma = molarcontantPressureSpecificHeat/molarcontantVolumeSpecificHeat;

    // calculate properties in cells attached to each boundary face

    vectorField momentum(faces_.size(), Zero);
    vectorField UCollected(faces_.size(), Zero);
    scalarField mass(faces_.size(), Zero);
    scalarField nParcels(faces_.size(), Zero);
    scalarField nParcelsInt(faces_.size(), Zero);
    vectorField velSqrMean(faces_.size(), Zero);
    vectorField velMeanSqr(faces_.size(), Zero);
    scalarField rotationalEnergy(faces_.size(), Zero);
    scalarField rotationalDof(faces_.size(), Zero);
    scalarField vDoF(faces_.size(), Zero);
    scalarField translationalTemperature(faces_.size(), Zero);
    scalarField rotationalTemperature(faces_.size(), Zero);
    scalarField vibrationalTemperature(faces_.size(), Zero);
    scalarField overallTemperature(faces_.size(), Zero);
    scalarField numberDensity(faces_.size(), Zero);
    scalarField massDensity(faces_.size(), Zero);
    scalarField pressure(faces_.size(), Zero);
    scalarField speedOfSound(faces_.size(), Zero);
    scalarField velocityCorrection(faces_.size(), Zero);
    scalarField massDensityCorrection(faces_.size(), Zero);

    const auto& cellOccupancy = cloud_.cellOccupancy();

    forAll(cells_, c)
    {
        const label cellI = cells_[c];
        const auto& parcelsInCell = cellOccupancy[cells_[c]];

        for (uspParcel* p : parcelsInCell)
        {
            label iD = typeIds_.find(p->typeId());

            if (iD != -1)
            {
                const scalar m =
                    cloud_.nParticle()
                   *cloud_.constProps(p->typeId()).mass();

                scalar CWF = cloud_.cellWF(cellI);
                scalar RWF = cloud_.axiRWF(p->position());

                momentum[c] += CWF*RWF*m*p->U();
                mass[c] += CWF*RWF*m;
            }

            nParcels[c] += 1.0;
            UCollected[c] += p->U();
            rotationalEnergy[c] += p->ERot();
            rotationalDof[c] += cloud_.constProps(p->typeId()).rotationalDoF();
            velSqrMean[c] += cmptSqr(p->U());
            velMeanSqr[c] += p->U();

            if (cloud_.constProps(p->typeId()).rotationalDoF() > VSMALL)
            {
                nParcelsInt[c] += 1.0;
            }
        }

        nTotalParcels_[c] += nParcels[c];

        nTotalParcelsInt_[c] += nParcelsInt[c];

        totalMomentum_[c] += momentum[c];

        totalMass_[c] += mass[c];

        UCollected_[c] += UCollected[c];

        totalRotationalEnergy_[c] += rotationalEnergy[c];

        totalRotationalDof_[c] += rotationalDof[c];

        massDensity[c] = totalMass_[c]/(cellVolume_[c]*nTimeSteps_);

        numberDensity[c] = massDensity[c]/molecularMass;

        if (nTotalParcels_[c] > 1)
        {
            velSqrMean_[c] += velSqrMean[c];
            velMeanSqr_[c] += velMeanSqr[c];

            totalVelSqrMean_[c] = velSqrMean_[c]/nTotalParcels_[c];
            totalVelMeanSqr_[c] = cmptSqr(velMeanSqr_[c]/nTotalParcels_[c]);

            translationalTemperature[c] =
                (0.5*molecularMass)*(2.0/(3.0*physicoChemical::k.value()))
               *(
                    cmptSum(totalVelSqrMean_[c])
                  - cmptSum(totalVelMeanSqr_[c])
                );

            UMean_[c] = UCollected_[c]/nTotalParcels_[c];

            if (translationalTemperature[c] < VSMALL)
            {
                translationalTemperature[c] = 300.00;
            }

            pressure[c] =
                numberDensity[c]
               *physicoChemical::k.value()
               *translationalTemperature[c];

            speedOfSound[c] =
                sqrt(gamma*gasConstant*translationalTemperature[c]);

            label faceI = faces_[c];
            vector n = mesh_.faceAreas()[faceI];
            n /= mag(n);

            if (nTimeSteps_ > 0)
            {
                massDensityCorrection[c] =
                    (outletPressure_ - pressure[c])
                   /(speedOfSound[c]*speedOfSound[c]);
                outletMassDensity_[c] =
                    massDensity[c] + massDensityCorrection[c];
            }
            else
            {
                outletMassDensity_[c] = massDensity[c];
            }

            // Liou and Fang, 2000, equation 26 STEP 1
            outletNumberDensity_[c] = outletMassDensity_[c]/molecularMass;

            outletTemperature_[c] =
                outletPressure_/(gasConstant*outletMassDensity_[c]);

            if (totalMass_[c] > 0)
            {
                outletVelocity_[c] = totalMomentum_[c]/totalMass_[c];
            }
            else
            {
                outletVelocity_[c] = vector::zero;
            }

            //velocity correction for each boundary cellI
            if (nTimeSteps_ > 0)
            {
                velocityCorrection[c] =
                    (pressure[c] - outletPressure_)
                   /(massDensity[c]*speedOfSound[c]);
                outletVelocity_[c] += velocityCorrection[c]*n;
            }

        }
    }
    
    // Compute number of parcels to insert
    computeParcelsToInsert
    (
        outletTemperature_,
        outletVelocity_,
        outletNumberDensity_,
        moleFractions_
    );
}


void Foam::uspLiouFangPressureOutlet::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uspLiouFangPressureOutlet::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uspGeneralBoundary::updateProperties(dict);

    setProperties();
}


void Foam::uspLiouFangPressureOutlet::setProperties()
{
    outletPressure_ = propsDict_.get<scalar>("outletPressure");

    // Read in the type ids
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // read in the mole fraction per specie

    const dictionary& moleFractionsDict(propsDict_.subDict("moleFractions"));

    moleFractions_.clear();

    moleFractions_.setSize(typeIds_.size(), 0.0);

    forAll(moleFractions_, i)
    {
        const word& moleculeName = cloud_.typeIdList()[typeIds_[i]];
        moleFractions_[i] = moleFractionsDict.get<scalar>(moleculeName);
    }

    // set the accumulator

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(), 0.0);
    }

    vibT_.setSize(typeIds_.size());

    forAll(vibT_, m)
    {
        vibT_[m].setSize(faces_.size(), 0.0);
    }

    vDof_.setSize(typeIds_.size());

    forAll(vDof_, m)
    {
        vDof_[m].setSize(faces_.size(), 0.0);
    }

    nTotalParcelsSpecies_.setSize(typeIds_.size());

    forAll(nTotalParcelsSpecies_, m)
    {
        nTotalParcelsSpecies_[m].setSize(faces_.size(), 0.0);
    }
}


// ************************************************************************* //
