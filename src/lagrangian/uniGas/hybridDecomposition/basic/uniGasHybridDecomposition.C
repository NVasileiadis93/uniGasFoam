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

#include "uniGasHybridDecomposition.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::uniGasHybridDecomposition::dictName("hybridDecompositionDict");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uniGasHybridDecomposition, 0);

    defineRunTimeSelectionTable(uniGasHybridDecomposition, dictionary);
};


Foam::uniGasHybridDecomposition::uniGasHybridDecomposition
(
    const Time& t,
    const polyMesh& mesh,
    uniGasCloud& owner
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
    decompositionInterval_(timeDict_.getOrDefault<label>("decompositionInterval",100)),
    resetAtDecomposition_(timeDict_.getOrDefault<bool>("resetAtDecomposition",true)),
    resetAtDecompositionUntilTime_(timeDict_.getOrDefault<scalar>("resetAtDecompositionUntilTime",VGREAT)),
    refinementPasses_(3),
    neighborLevels_(3),
    maxNeighborFraction_(0.4),
    boundCoeff_(0.5),
    neighborCells_()
{

}

Foam::autoPtr<Foam::uniGasHybridDecomposition> Foam::uniGasHybridDecomposition::New
(
    const Time& t,
    const polyMesh& mesh,
    uniGasCloud& owner
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

    Info<< "Selecting hybridDecompositionModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "uniGasHybridDecomposition",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<uniGasHybridDecomposition>(cstrIter()(t, mesh, owner));
}

// ************************************************************************* //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasHybridDecomposition::fetchCellNeighborhood
(
    const label cell,
    const label nLevels,
    DynamicList<label>& neighborCells
)
{

    // initialize neighborCells
    neighborCells.setCapacity(0);

    // get cells in extended neighbor area
    label neighborCandidate;
    label neighborCellinitial = 0;
    label neighborCellFinal = 0;
    label neighborCellCapacity = 1;
    boolList isNeighbor(mesh_.nCells());
    isNeighbor = false;
    isNeighbor[cell]=true;
    neighborCells.append(cell);

    for(label level=1; level<=nLevels; level++)
    {

        neighborCellinitial=neighborCellFinal;
        neighborCellFinal=neighborCellCapacity-1;

        for(label cellI=neighborCellinitial; cellI<=neighborCellFinal; cellI++)
        {

            forAll(mesh_.cellCells()[neighborCells[cellI]], cellJ)
            {

                  neighborCandidate=mesh_.cellCells()[neighborCells[cellI]][cellJ];

                  if (!isNeighbor[neighborCandidate])
                  {

                      neighborCellCapacity++;
                      isNeighbor[neighborCandidate]=true;
                      neighborCells.append(neighborCandidate);

                  }
            }
        }
    }

}

void Foam::uniGasHybridDecomposition::updateProperties()
{

    timeDict_ = hybridDecompositionDict_.subDict("timeProperties");

    timeDict_.readIfPresent("decompositionInterval", decompositionInterval_);

    timeDict_.readIfPresent("resetAtDecomposition", resetAtDecomposition_);

    timeDict_.readIfPresent("resetAtDecompositionUntilTime", resetAtDecompositionUntilTime_);

}

// ************************************************************************* //
