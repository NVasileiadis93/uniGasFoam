/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "collisionPartnerSelection.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(collisionPartnerSelection, 0);
defineRunTimeSelectionTable(collisionPartnerSelection, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::collisionPartnerSelection::collisionPartnerSelection
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    rndGen_(cloud_.rndGen())
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::collisionPartnerSelection>
Foam::collisionPartnerSelection::New
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("collisionPartnerSelectionModel"));

    Info<< "Selecting collisionPartnerSelectionModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "collisionPartnerSelectionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<collisionPartnerSelection>(cstrIter()(mesh, cloud, dict));
}


// ************************************************************************* //
