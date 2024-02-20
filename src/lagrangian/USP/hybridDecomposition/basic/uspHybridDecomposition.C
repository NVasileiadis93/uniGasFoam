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


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::uspHybridDecomposition::dictName("hybridDecompositionDict");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uspHybridDecomposition, 0);

    defineRunTimeSelectionTable(uspHybridDecomposition, dictionary);
};


Foam::uspHybridDecomposition::uspHybridDecomposition
(
    const Time& t,
    const polyMesh& mesh,
    uspCloud& owner
)
:
    time_(t),
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(owner),
    hybridDecompositionDict_
    (
        IOobject
        (
            dictName,
            t.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    timeDict_(hybridDecompositionDict_.subDict("timeProperties")),
    decompositionInterval_(timeDict_.get<label>("decompositionInterval")),
    resetAtDecomposition_(timeDict_.getOrDefault<bool>("resetAtDecomposition",true)),
    resetAtDecompositionUntilTime_(timeDict_.getOrDefault<scalar>("resetAtDecompositionUntilTime",VGREAT)),
    refinementPasses_(10),
    boundCoeff_(0.5)
{

}

Foam::autoPtr<Foam::uspHybridDecomposition> Foam::uspHybridDecomposition::New
(
    const Time& t,
    const polyMesh& mesh,
    uspCloud& owner
)
{

    IOdictionary dict
    (
        IOobject
        (
            dictName,
            t.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

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

    return autoPtr<uspHybridDecomposition>(cstrIter()(t, mesh, owner));
}

// ************************************************************************* //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspHybridDecomposition::updateProperties()
{

    timeDict_ = hybridDecompositionDict_.subDict("timeProperties");

    timeDict_.readIfPresent("decompositionInterval", decompositionInterval_);

    timeDict_.readIfPresent("resetAtDecomposition", resetAtDecomposition_);

    timeDict_.readIfPresent("resetAtDecompositionUntilTime", resetAtDecompositionUntilTime_);

}

// ************************************************************************* //
