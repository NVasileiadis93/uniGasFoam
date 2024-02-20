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

#include "uspHybridDecomposition.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uspHybridDecomposition, 0);

    defineRunTimeSelectionTable(uspHybridDecomposition, dictionary);
};


Foam::uspHybridDecomposition::uspHybridDecomposition
(
    const dictionary& dict,
    const polyMesh& mesh,
    uspCloud& owner
)
:
    dict_(dict),
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(owner),
    decompositionInterval_(dict.subDict("decompositionProperties").get<label>("decompositionInterval")),
    breakdownMax_(dict.subDict("decompositionProperties").get<scalar>("breakdownMax")),
    theta_(dict.subDict("decompositionProperties").getOrDefault<scalar>("theta",1.0)),
    smoothingPasses_(dict.subDict("decompositionProperties").getOrDefault<scalar>("smoothingPasses",0)),
    refinementPasses_(10),
    boundCoeff_(0.5)
{

}

Foam::autoPtr<Foam::uspHybridDecomposition> Foam::uspHybridDecomposition::New
(
    const dictionary& dict,
    const polyMesh& mesh,
    uspCloud& owner
)
{
    const word modelType(dict.get<word>("decompositionModel"));

    Info<< "Selecting hybrid decomposition model " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "uspHybridDecomposition",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<uspHybridDecomposition>(cstrIter()(dict, mesh, owner));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::uspHybridDecomposition::dict() const
{
    return dict_;
}


// ************************************************************************* //
