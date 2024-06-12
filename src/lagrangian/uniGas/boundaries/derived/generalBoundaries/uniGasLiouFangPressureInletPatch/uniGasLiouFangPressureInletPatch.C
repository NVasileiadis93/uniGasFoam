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

\*---------------------------------------------------------------------------*/

#include "uniGasLiouFangPressureInletPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{
defineTypeNameAndDebug(uniGasLiouFangPressureInletPatch, 0);

addToRunTimeSelectionTable
(
    uniGasGeneralBoundary,
    uniGasLiouFangPressureInletPatch,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uniGasLiouFangPressureInletPatch::uniGasLiouFangPressureInletPatch
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    inletNumberDensity_(),
    moleFractions_(),
    inletPressure_(propsDict_.get<scalar>("inletPressure")),
    inletTemperature_(propsDict_.get<scalar>("inletTemperature")),
    inletVelocity_(faces_.size(), Zero),
    theta_(propsDict_.getOrDefault<scalar>("theta",1.0)),
    previousInletVelocity_(faces_.size(), Zero)
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    if (0.0 > theta_ || theta_ > 1.0)
    {
        FatalErrorInFunction
            << "Theta must be a value between 0 and 1 " << nl << "in: "
            << mesh_.time().system()/uniGasBoundaries::dictName
            << exit(FatalError);
    }

    // Get type IDs
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // Set the accumulator
    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(), 0.0);
    }

    // Read in the mole fraction per specie
    const dictionary& moleFractionsDict(propsDict_.subDict("moleFractions"));

    moleFractions_.clear();

    moleFractions_.setSize(typeIds_.size(), 0.0);

    forAll(moleFractions_, i)
    {
        const word& moleculeName = cloud_.typeIdList()[typeIds_[i]];
        moleFractions_[i] = moleFractionsDict.get<scalar>(moleculeName);
    }

    inletNumberDensity_ = inletPressure_ / (physicoChemical::k.value()*inletTemperature_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uniGasLiouFangPressureInletPatch::initialConfiguration()
{}


void uniGasLiouFangPressureInletPatch::calculateProperties()
{}


void uniGasLiouFangPressureInletPatch::controlParcelsBeforeMove()
{
    insertParcels
    (
        inletTemperature_,
        inletVelocity_
    );

    previousInletVelocity_ = inletVelocity_;
}


void uniGasLiouFangPressureInletPatch::controlParcelsBeforeCollisions()
{}


void uniGasLiouFangPressureInletPatch::controlParcelsAfterCollisions()
{

    vectorField momentum(faces_.size(), Zero);
    vectorField newInletVelocity(faces_.size(), Zero);
    scalarField mass(faces_.size(), scalar(0));

    const List<DynamicList<uniGasParcel*>>& cellOccupancy = cloud_.cellOccupancy();

    forAll(cells_, c)
    {
        const label cellI = cells_[c];
        const List<uniGasParcel*>& parcelsInCell = cellOccupancy[cells_[c]];

        forAll(parcelsInCell, pIC)
        {

            uniGasParcel* p = parcelsInCell[pIC];

            const scalar m = cloud_.nParticle()*cloud_.constProps(p->typeId()).mass();

            scalar CWF = cloud_.cellWF(cellI);
            scalar RWF = cloud_.axiRWF(p->position());

            momentum[c] += CWF*RWF*m*p->U();
            mass[c] += CWF*RWF*m;

        }

        if (mass[c]>0)
        {
            newInletVelocity[c] = momentum[c]/mass[c];
        }
        
        inletVelocity_[c] = theta_*newInletVelocity[c] + (1.0 - theta_)*previousInletVelocity_[c];
    }

    // Compute number of parcels to insert
    computeParcelsToInsert
    (
        inletNumberDensity_,
        moleFractions_,
        inletTemperature_,
        inletVelocity_
    );
}

void uniGasLiouFangPressureInletPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void uniGasLiouFangPressureInletPatch::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uniGasGeneralBoundary::updateProperties(dict);
}

} // End namespace Foam

// ************************************************************************* //
