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

#include "dsmcCollisionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dsmcCollisionModel, 0);

    defineRunTimeSelectionTable(dsmcCollisionModel, dictionary);
};


Foam::dsmcCollisionModel::dsmcCollisionModel(uniGasCloud& owner)
:
    dict_(dictionary::null),
    cloud_(owner)
{}


Foam::dsmcCollisionModel::dsmcCollisionModel
(
    const dictionary& dict,
    uniGasCloud& owner
)
:
    dict_(dict),
    cloud_(owner)
{}


Foam::autoPtr<Foam::dsmcCollisionModel> Foam::dsmcCollisionModel::New
(
    const dictionary& dict,
    uniGasCloud& owner
)
{
    const word modelType(dict.get<word>("dsmcCollisionModel"));

    Info<< "Selecting dsmcCollisionModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "dsmcCollisionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<dsmcCollisionModel>(cstrIter()(dict, owner));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::dsmcCollisionModel::dict() const
{
    return dict_;
}


// ************************************************************************* //
