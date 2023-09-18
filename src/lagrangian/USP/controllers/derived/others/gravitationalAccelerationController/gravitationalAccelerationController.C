/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "gravitationalAccelerationController.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(gravitationalAccelerationController, 0);

addToRunTimeSelectionTable
(
    uspStateController,
    gravitationalAccelerationController,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gravitationalAccelerationController::gravitationalAccelerationController
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspStateController(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    acceleration_(propsDict_.get<vector>("acceleration"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    singleValueController() = true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gravitationalAccelerationController::initialConfiguration()
{}


void Foam::gravitationalAccelerationController::calculateProperties()
{}


void Foam::gravitationalAccelerationController::controlParcelsBeforeMove()
{
    forAll(controlZone(), c)
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();
        const label celli = controlZone()[c];

        for (uspParcel* p : cellOccupancy[celli])
        {
            p->U() += 0.5*acceleration_*mesh_.time().deltaTValue();
        }
    }
}


void Foam::gravitationalAccelerationController::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    const Time& runTime = mesh_.time();

    if (runTime.writeTime())
    {
        NotImplemented;
    }
}


void Foam::gravitationalAccelerationController::controlParcelsBeforeCollisions()
{
    forAll(controlZone(), c)
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();
        const label celli = controlZone()[c];

        for (uspParcel* p : cellOccupancy[celli])
        {
            p->U() -= 0.5*acceleration_*mesh_.time().deltaTValue();
        }
    }
}


void Foam::gravitationalAccelerationController::controlParcelsAfterCollisions()
{
    forAll(controlZone(), c)
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();
        const label celli = controlZone()[c];

        for (uspParcel* p : cellOccupancy[celli])
        {
            p->U() += acceleration_*mesh_.time().deltaTValue();
        }
    }
}


void Foam::gravitationalAccelerationController::updateProperties
(
    const dictionary& dict
)
{
    uspStateController::updateProperties(dict);
    propsDict_ = dict.subDict(typeName + "Properties");
    setProperties();
}


void Foam::gravitationalAccelerationController::setProperties()
{
    acceleration_ = propsDict_.get<vector>("acceleration");
}


// ************************************************************************* //
