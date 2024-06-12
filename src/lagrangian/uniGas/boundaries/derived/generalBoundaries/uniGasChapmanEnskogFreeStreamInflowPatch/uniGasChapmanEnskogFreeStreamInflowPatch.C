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

#include "uniGasChapmanEnskogFreeStreamInflowPatch.H"
#include "mathematicalConstants.H"
#include "Random.H"
#include "uniGasCloud.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{
defineTypeNameAndDebug(uniGasChapmanEnskogFreeStreamInflowPatch, 0);

addToRunTimeSelectionTable
(
    uniGasGeneralBoundary,
    uniGasChapmanEnskogFreeStreamInflowPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasChapmanEnskogFreeStreamInflowPatch::uniGasChapmanEnskogFreeStreamInflowPatch
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    translationalTemperature_(propsDict_.get<scalar>("translationalTemperature")),
    rotationalTemperature_(propsDict_.get<scalar>("rotationalTemperature")),
    vibrationalTemperature_(propsDict_.get<scalar>("vibrationalTemperature")),
    electronicTemperature_(propsDict_.get<scalar>("electronicTemperature")),
    velocity_(propsDict_.get<vector>("velocity")),
    heatFlux_(propsDict_.get<vector>("heatFlux")),
    stress_(propsDict_.get<tensor>("stress"))
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    // Get type IDs
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // Set the accumulator
    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(), Zero);
    }

    // Read in the number density per specie
    const dictionary& numberDensitiesDict(propsDict_.subDict("numberDensities"));

    numberDensities_.clear();

    numberDensities_.setSize(typeIds_.size(), Zero);

    forAll(numberDensities_, i)
    {
        const word& moleculeName = cloud_.typeIdList()[typeIds_[i]];
        numberDensities_[i] = numberDensitiesDict.get<scalar>(moleculeName);
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasChapmanEnskogFreeStreamInflowPatch::initialConfiguration()
{}


void Foam::uniGasChapmanEnskogFreeStreamInflowPatch::calculateProperties()
{}


void Foam::uniGasChapmanEnskogFreeStreamInflowPatch::controlParcelsBeforeMove()
{
    computeParcelsToInsert
    (
        numberDensities_,
        translationalTemperature_,
        velocity_,
        heatFlux_,
        stress_
    );

    insertParcels
    (
        numberDensities_,
        translationalTemperature_,
        rotationalTemperature_,
        vibrationalTemperature_,
        electronicTemperature_,
        velocity_,
        heatFlux_,
        stress_
    );
}


void Foam::uniGasChapmanEnskogFreeStreamInflowPatch::controlParcelsBeforeCollisions()
{}


void Foam::uniGasChapmanEnskogFreeStreamInflowPatch::controlParcelsAfterCollisions()
{}


void Foam::uniGasChapmanEnskogFreeStreamInflowPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uniGasChapmanEnskogFreeStreamInflowPatch::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uniGasGeneralBoundary::updateProperties(dict);
}

// ************************************************************************* //
