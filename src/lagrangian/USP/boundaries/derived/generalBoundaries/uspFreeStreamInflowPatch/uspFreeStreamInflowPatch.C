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

#include "uspFreeStreamInflowPatch.H"
#include "mathematicalConstants.H"
#include "Random.H"
#include "uspCloud.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{
defineTypeNameAndDebug(uspFreeStreamInflowPatch, 0);

addToRunTimeSelectionTable
(
    uspGeneralBoundary,
    uspFreeStreamInflowPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspFreeStreamInflowPatch::uspFreeStreamInflowPatch
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    translationalTemperature_(propsDict_.get<scalar>("translationalTemperature")),
    rotationalTemperature_(propsDict_.get<scalar>("rotationalTemperature")),
    vibrationalTemperature_(propsDict_.get<scalar>("vibrationalTemperature")),
    electronicTemperature_(propsDict_.get<scalar>("electronicTemperature")),
    velocity_(propsDict_.get<vector>("velocity"))
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
    const dictionary& numberDensitiesDict
    (
        propsDict_.subDict("numberDensities")
    );

    numberDensities_.clear();

    numberDensities_.setSize(typeIds_.size(), Zero);

    forAll(numberDensities_, i)
    {
        const word& moleculeName = cloud_.typeIdList()[typeIds_[i]];
        numberDensities_[i] = numberDensitiesDict.get<scalar>(moleculeName);
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspFreeStreamInflowPatch::initialConfiguration()
{}


void Foam::uspFreeStreamInflowPatch::calculateProperties()
{}


void Foam::uspFreeStreamInflowPatch::controlParcelsBeforeMove()
{
    computeParcelsToInsert
    (
        numberDensities_,
        translationalTemperature_,
        velocity_
    );

    insertParcels
    (
        translationalTemperature_,
        rotationalTemperature_,
        vibrationalTemperature_,
        electronicTemperature_,
        velocity_
    );
}


void Foam::uspFreeStreamInflowPatch::controlParcelsBeforeCollisions()
{}


void Foam::uspFreeStreamInflowPatch::controlParcelsAfterCollisions()
{}


void Foam::uspFreeStreamInflowPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uspFreeStreamInflowPatch::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uspGeneralBoundary::updateProperties(dict);
}

// ************************************************************************* //
