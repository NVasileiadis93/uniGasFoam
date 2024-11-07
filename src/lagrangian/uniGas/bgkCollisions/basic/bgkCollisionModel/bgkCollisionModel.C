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

#include "bgkCollisionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bgkCollisionModel, 0);

    defineRunTimeSelectionTable(bgkCollisionModel, dictionary);
};


Foam::bgkCollisionModel::bgkCollisionModel
(
    const dictionary& dict,
    const polyMesh& mesh,
    uniGasCloud& owner
)
:
    dict_(dict),
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(owner)
{

}

Foam::autoPtr<Foam::bgkCollisionModel> Foam::bgkCollisionModel::New
(
    const dictionary& dict,
    const polyMesh& mesh,
    uniGasCloud& owner
)
{
    const word modelType(dict.get<word>("bgkCollisionModel"));

    Info<< "Selecting bgkCollisionModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "bgkCollisionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<bgkCollisionModel>(cstrIter()(dict, mesh, owner));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::bgkCollisionModel::dict() const
{
    return dict_;
}


// ************************************************************************* //
