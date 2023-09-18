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

#include "relaxationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relaxationModel, 0);

    defineRunTimeSelectionTable(relaxationModel, dictionary);
};


Foam::relaxationModel::relaxationModel
(
    const dictionary& dict,
    const polyMesh& mesh,
    uspCloud& owner
)
:
    dict_(dict),
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(owner)
{

}

Foam::autoPtr<Foam::relaxationModel> Foam::relaxationModel::New
(
    const dictionary& dict,
    const polyMesh& mesh,
    uspCloud& owner
)
{
    const word modelType(dict.get<word>("relaxationModel"));

    Info<< "Selecting relaxationModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "relaxationModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<relaxationModel>(cstrIter()(dict, mesh, owner));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::relaxationModel::dict() const
{
    return dict_;
}


// ************************************************************************* //
