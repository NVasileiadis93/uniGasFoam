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

#include "uspCyclicBoundary.H"
#include "processorCyclicPolyPatch.H"
#include "uspBoundaries.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspCyclicBoundary, 0);
defineRunTimeSelectionTable(uspCyclicBoundary, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspCyclicBoundary::uspCyclicBoundary
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspBoundaryBase(mesh, cloud, dict, "cyclic"),
    neighbPatchName_(),
    neighbPatchId_(0),

    faces_(),
    coupledFacesA_(),
    coupledFacesB_(),
    cellsA_(),
    cellsB_()
{
    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

    if (!isA<cyclicPolyPatch>(patch))
    {
        FatalErrorInFunction
            << "Patch: " << patchName_ << " is not a cyclic boundary. "
            << nl << "in: "
            << mesh_.time().system()/uspBoundaries::dictName
            << exit(FatalError);
    }

    const cyclicPolyPatch& cyclicPatch = refCast<const cyclicPolyPatch>(patch);

    neighbPatchName_ = cyclicPatch.neighbPatchName();

    neighbPatchId_ = cyclicPatch.neighbPatchID();

    getCoupledFaces(mesh_);
}


void Foam::uspCyclicBoundary::getCoupledFaces(const polyMesh& mesh)
{
    // In parallel, openFOAM replaces cyclic boundaries with processor
    // boundaries if processor boundaries coincide with cyclic boundaries.
    // We therefore need to do the following corrections to get the
    // correct faces on the selected cyclic patch.

    if (Pstream::parRun())
    {
        DynamicList<label> coupledFacesA(0);
        DynamicList<label> coupledFacesB(0);
        DynamicList<label> cellsA(0);
        DynamicList<label> cellsB(0);

        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

            for (label i = 0; i < patch.size(); ++i)
            {
                label globalFaceI = patch.start() + i;
                label cellI = patch.faceCells()[i];
                coupledFacesA.append(globalFaceI);
                cellsA.append(cellI);
            }

            const polyPatch& patchN = mesh_.boundaryMesh()[neighbPatchId_];

            for (label i = 0; i < patchN.size(); ++i)
            {
                label globalFaceI = patchN.start() + i;
                label cellI = patchN.faceCells()[i];
                coupledFacesB.append(globalFaceI);
                cellsB.append(cellI);
            }
        }

        // Check for processor-cyclic poly patches

        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];

            if (isA<processorCyclicPolyPatch>(patch))
            {
                const word& patchName =
                    refCast<const processorCyclicPolyPatch>
                    (
                        patch
                    ).referPatchName();

                if (patchName == patchName_)
                {
                    for (label i = 0; i < patch.size(); ++i)
                    {
                        label globalFaceI = patch.start() + i;
                        label cellI = patch.faceCells()[i];
                        coupledFacesA.append(globalFaceI);
                        cellsA.append(cellI);
                    }
                }
                else if (patchName == neighbPatchName_)
                {
                    for (label i = 0; i < patch.size(); ++i)
                    {
                        label globalFaceI = patch.start() + i;
                        label cellI = patch.faceCells()[i];
                        coupledFacesB.append(globalFaceI);
                        cellsB.append(cellI);
                    }
                }
            }
        }

        coupledFacesA.shrink();
        coupledFacesB.shrink();
        cellsA.shrink();
        cellsB.shrink();

        coupledFacesA_.transfer(coupledFacesA);
        coupledFacesB_.transfer(coupledFacesB);
        cellsA_.transfer(cellsA);
        cellsB_.transfer(cellsB);
    }
    else
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

        coupledFacesA_.setSize(patch.size(), -1);
        cellsA_.setSize(patch.size(), -1);

        for (label i = 0; i < patch.size(); ++i)
        {
            label globalFaceI = patch.start() + i;
            label cellI = patch.faceCells()[i];
            coupledFacesA_[i] = globalFaceI;
            cellsA_[i] = cellI;
        }

        const polyPatch& patchN = mesh_.boundaryMesh()[neighbPatchId_];

        coupledFacesB_.setSize(patchN.size(), -1);
        cellsB_.setSize(patchN.size(), -1);

        for (label i = 0; i < patchN.size(); ++i)
        {
            label globalFaceI = patchN.start() + i;
            label cellI = patchN.faceCells()[i];
            coupledFacesB_[i] = globalFaceI;
            cellsB_[i] = cellI;
        }
    }

    label totalFaces = coupledFacesA_.size() + coupledFacesB_.size();

    faces_.setSize(totalFaces, -1);

    label index = 0;

    forAll(coupledFacesA_, fA)
    {
        faces_[index] = coupledFacesA_[fA];
        ++index;
    }

    forAll(coupledFacesB_, fB)
    {
        faces_[index] = coupledFacesB_[fB];
        ++index;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::uspCyclicBoundary> Foam::uspCyclicBoundary::New
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
{
    word modelType(dict.get<word>("boundaryModel"));

    Info<< "Selecting uspCyclicBoundaryModel " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "uspCyclicBoundary",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<uspCyclicBoundary>(cstrIter()(mesh, cloud, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::uspCyclicBoundary::controlPatch() const
{
    return coupledFacesA_;
}


const Foam::labelList& Foam::uspCyclicBoundary::controlZone() const
{
    return cellsA_;
}


const Foam::labelList& Foam::uspCyclicBoundary::allFaces() const
{
    return faces_;
}


// ************************************************************************* //
