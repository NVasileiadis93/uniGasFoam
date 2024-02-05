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

#include "binaryCollisionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(binaryCollisionModel, 0);

    defineRunTimeSelectionTable(binaryCollisionModel, dictionary);
};


Foam::binaryCollisionModel::binaryCollisionModel(uspCloud& owner)
:
    dict_(dictionary::null),
    cloud_(owner)
{}


Foam::binaryCollisionModel::binaryCollisionModel
(
    const dictionary& dict,
    uspCloud& owner
)
:
    dict_(dict),
    cloud_(owner)
{}


Foam::autoPtr<Foam::binaryCollisionModel> Foam::binaryCollisionModel::New
(
    const dictionary& dict,
    uspCloud& owner
)
{
    const word modelType(dict.get<word>("binaryCollisionModel"));

    Info<< "Selecting binaryCollisionModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "binaryCollisionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<binaryCollisionModel>(cstrIter()(dict, owner));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::binaryCollisionModel::dict() const
{
    return dict_;
}


// ************************************************************************* //
