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

#include "noDSMCCollision.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(noDSMCCollision, 0);
    addToRunTimeSelectionTable
    (
        dsmcCollisionModel,
        noDSMCCollision,
        dictionary
    );
};


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::noDSMCCollision::noDSMCCollision
(
    const dictionary& dict,
    uniGasCloud& cloud
)
:
    dsmcCollisionModel(cloud)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::noDSMCCollision::active() const
{
    return false;
}


Foam::scalar Foam::noDSMCCollision::sigmaTcR
(
    const uniGasParcel& pP,
    const uniGasParcel& pQ
) const
{
    FatalErrorInFunction
        << "sigmaTcR called on noDSMCCollision model, this should "
        << "not happen, this is not an actual model." << nl
        << "Enclose calls to sigmaTcR within a if (dsmcCollision().active()) "
        << " check."
        << abort(FatalError);

    return 0.0;
}


void Foam::noDSMCCollision::collide
(
    uniGasParcel& pP,
    uniGasParcel& pQ,
    label& cellI
)
{}


const Foam::dictionary& Foam::noDSMCCollision::propertiesDict() const
{
    return propertiesDict_;
}


// ************************************************************************* //
