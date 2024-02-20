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

#include "uspBoundaries.H"
#include "emptyPolyPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::uspBoundaries::dictName("boundariesDict");


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

uspBoundaries::uspBoundaries(const Time& t, const polyMesh& mesh)
:
    time_(t),
    uspBoundariesDict_
    (
        IOobject
        (
            uspBoundaries::dictName,
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),

    nPatchBoundaryModels_(0),
    nCyclicBoundaryModels_(0),
    nGeneralBoundaryModels_(0),

    patchBoundaryList_(),
    patchBoundaryNames_(),
    patchBoundaryIds_(),
    pBFixedPathNames_(),
    patchBoundaryModels_(),
    patchToModelId_(mesh.boundaryMesh().size(), -1),

    cyclicBoundaryList_(),
    cyclicBoundaryNames_(),
    cyclicBoundaryIds_(),
    cMFixedPathNames_(),
    cyclicBoundaryModels_(),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1),

    generalBoundaryList_(),
    generalBoundaryNames_(),
    generalBoundaryIds_(),
    gMFixedPathNames_(),
    generalBoundaryModels_()
{}


uspBoundaries::uspBoundaries
(
    const Time& t,
    const polyMesh& mesh,
    uspCloud& cloud
)
:
    time_(t),
    uspBoundariesDict_
    (
        IOobject
        (
            uspBoundaries::dictName,
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    nPatchBoundaryModels_(0),
    nCyclicBoundaryModels_(0),
    nGeneralBoundaryModels_(0),

    patchBoundaryList_(uspBoundariesDict_.lookup("uspPatchBoundaries")),
    patchBoundaryNames_(patchBoundaryList_.size()),
    patchBoundaryIds_(patchBoundaryList_.size()),
    pBFixedPathNames_(patchBoundaryList_.size()),
    patchBoundaryModels_(patchBoundaryList_.size()),
    patchToModelId_(mesh.boundaryMesh().size(), -1),

    cyclicBoundaryList_(uspBoundariesDict_.lookup("uspCyclicBoundaries")),
    cyclicBoundaryNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryIds_(cyclicBoundaryList_.size()),
    cMFixedPathNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryModels_(cyclicBoundaryList_.size()),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1),

    generalBoundaryList_(uspBoundariesDict_.lookup("uspGeneralBoundaries")),
    generalBoundaryNames_(generalBoundaryList_.size()),
    generalBoundaryIds_(generalBoundaryList_.size()),
    gMFixedPathNames_(generalBoundaryList_.size()),
    generalBoundaryModels_(generalBoundaryList_.size())
{
    Info << "Creating the boundary models: " << nl << endl;

    // Patch boundaries

    if (patchBoundaryModels_.size() > 0)
    {
        forAll(patchBoundaryModels_, p)
        {
            const entry& boundaryI = patchBoundaryList_[p];
            const dictionary& boundaryIDict = boundaryI.dict();

            patchBoundaryModels_[p] = autoPtr<uspPatchBoundary>
            (
                uspPatchBoundary::New(mesh, cloud, boundaryIDict)
            );

            patchBoundaryNames_[p] = patchBoundaryModels_[p]->type();
            patchBoundaryIds_[p] = p;
            ++nPatchBoundaryModels_;
        }
    }

    checkPatchBoundaryModels(mesh);

    // Cyclic boundaries

    if (cyclicBoundaryModels_.size() > 0)
    {
        forAll(cyclicBoundaryModels_, c)
        {
            const entry& boundaryI = cyclicBoundaryList_[c];
            const dictionary& boundaryIDict = boundaryI.dict();

            cyclicBoundaryModels_[c] = autoPtr<uspCyclicBoundary>
            (
                uspCyclicBoundary::New(mesh, cloud, boundaryIDict)
            );

            cyclicBoundaryNames_[c] = cyclicBoundaryModels_[c]->type();
            cyclicBoundaryIds_[c] = c;
            ++nCyclicBoundaryModels_;
        }
    }

    checkCyclicBoundaryModels(mesh);

    // General boundaries

    if (generalBoundaryModels_.size() > 0)
    {
        forAll(generalBoundaryModels_, g)
        {
            const entry& boundaryI = generalBoundaryList_[g];
            const dictionary& boundaryIDict = boundaryI.dict();

            generalBoundaryModels_[g] = autoPtr<uspGeneralBoundary>
            (
                uspGeneralBoundary::New(mesh, cloud, boundaryIDict)
            );

            generalBoundaryNames_[g] = generalBoundaryModels_[g]->type();
            generalBoundaryIds_[g] = g;
            ++nGeneralBoundaryModels_;
        }
    }

    // Creating directories
    if (nPatchBoundaryModels_ > 0)
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if (!isDir(boundariesPath))
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/usp
        fileName uspBoundariesPath(boundariesPath/"usp");

        if (!isDir(uspBoundariesPath))
        {
            mkDir(uspBoundariesPath);
        }

        // directory: case/boundaries/usp/patchBoundaryModels
        fileName patchBoundaryModelsPath
        (
            uspBoundariesPath/"patchBoundaryModels"
        );

        if (!isDir(patchBoundaryModelsPath))
        {
            mkDir(patchBoundaryModelsPath);
        }

        forAll(patchBoundaryModels_, p)
        {
            if (patchBoundaryModels_[p]->writeInCase())
            {
                // directory: case/boundaries/usp/
                // patchBoundaryModels/<patchBoundaryModel>
                fileName patchBoundaryModelPath
                (
                    patchBoundaryModelsPath/patchBoundaryNames_[p]
                );

                if (!isDir(patchBoundaryModelPath))
                {
                    mkDir(patchBoundaryModelPath);
                }

                const word& patchName = patchBoundaryModels_[p]->patchName();

                // directory: case/boundaries/usp/
                // patchBoundaryModels/<patchBoundaryModel>/<patchName>

                fileName patchPath(patchBoundaryModelPath/patchName);

                if (!isDir(patchPath))
                {
                    mkDir(patchPath);
                }

                pBFixedPathNames_[p] = patchPath;
            }
        }
    }


    // Creating directories
    if (nCyclicBoundaryModels_ > 0)
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if (!isDir(boundariesPath))
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/usp
        fileName uspBoundariesPath(boundariesPath/"usp");

        if (!isDir(uspBoundariesPath))
        {
            mkDir(uspBoundariesPath);
        }

        // directory: case/boundaries/usp/cyclicBoundaryModels
        fileName cyclicBoundaryModelsPath
                    (uspBoundariesPath/"cyclicBoundaryModels");

        if (!isDir(cyclicBoundaryModelsPath))
        {
            mkDir(cyclicBoundaryModelsPath);
        }

        forAll(cyclicBoundaryModels_, c)
        {
            if (cyclicBoundaryModels_[c]->writeInCase())
            {
                // directory: case/boundaries/usp/
                // cyclicBoundaryModels/<cyclicBoundaryModel>
                fileName cyclicBoundaryModelPath
                (
                    cyclicBoundaryModelsPath/cyclicBoundaryNames_[c]
                );

                if (!isDir(cyclicBoundaryModelPath))
                {
                    mkDir(cyclicBoundaryModelPath);
                }

                const word& patchName = cyclicBoundaryModels_[c]->patchName();

                // directory: case/boundaries/usp
                // cyclicBoundaryModels/<cyclicBoundaryModel>/<patchName>

                fileName patchPath(cyclicBoundaryModelPath/patchName);

                if (!isDir(patchPath))
                {
                    mkDir(patchPath);
                }

                cMFixedPathNames_[c] = patchPath;
            }
        }
    }

    // Creating directories
    if (nGeneralBoundaryModels_ > 0)
    {
        // directory: case/boundaries
        fileName boundariesPath(time_.path()/"boundaries");

        if (!isDir(boundariesPath))
        {
            mkDir(boundariesPath);
        }

        // directory: case/boundaries/usp
        fileName uspBoundariesPath(boundariesPath/"usp");

        if (!isDir(uspBoundariesPath))
        {
            mkDir(uspBoundariesPath);
        }

        // directory: case/boundaries/usp/cyclicBoundaryModels
        fileName generalBoundaryModelsPath
        (
            uspBoundariesPath/"generalBoundaryModels"
        );

        if (!isDir(generalBoundaryModelsPath))
        {
            mkDir(generalBoundaryModelsPath);
        }

        forAll(generalBoundaryModels_, g)
        {
            if (generalBoundaryModels_[g]->writeInCase())
            {
                // directory: case/boundaries/usp
                // generalBoundaryModels/<generalBoundaryModel>
                fileName generalBoundaryModelPath
                (
                    generalBoundaryModelsPath/generalBoundaryNames_[g]
                );

                if (!isDir(generalBoundaryModelPath))
                {
                    mkDir(generalBoundaryModelPath);
                }

                const word& patchName = generalBoundaryModels_[g]->patchName();

                // directory: case/boundaries/usp
                // generalBoundaryModels/<generalBoundaryModel>/<patchName>

                fileName patchPath(generalBoundaryModelPath/patchName);

                if (!isDir(patchPath))
                {
                    mkDir(patchPath);
                }

                gMFixedPathNames_[g] = patchPath;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uspBoundaries::checkCyclicBoundaryModels(const polyMesh& mesh)
{
    label nPolyPatches = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (isA<cyclicPolyPatch>(patch))
        {
            label patchIndex = patch.index();

            forAll(cyclicBoundaryModels_, c)
            {
                const label patchId = cyclicBoundaryModels_[c]->patchId();

                if (patchIndex == patchId)
                {
                    ++nPolyPatches;
                    cyclicBoundaryToModelId_[patchi] = c;
                }
            }
        }
    }

    if (nPolyPatches != nCyclicBoundaryModels_)
    {
        FatalErrorInFunction
            << " Number of cyclic boundary models = " << nCyclicBoundaryModels_
            << " chosen in the boundaryiesDict are inconsistent."
            << abort(FatalError);
    }
}


void uspBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)
{
    // Check that all poly-patches defined within blockMeshDict,
    // each have one model.

    label nPolyPatches = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (!polyPatch::constraintType(patch.type()))
        {
            ++nPolyPatches;

            label patchIndex = patch.index();

            label nPatches = 0;

            forAll(patchBoundaryModels_, p)
            {
                const label patchId = patchBoundaryModels_[p]->patchId();

                if (patchIndex == patchId)
                {
                    ++nPatches;
                    patchToModelId_[patchi] = p;
                }
            }

            if (nPatches > 1)
            {
                FatalIOErrorInFunction(uspBoundariesDict_)
                    << " Only one patch boundary model per poly-patch, [name: "
                    << patch.name()
                    << "]. Number of models chosen for this patch are: "
                    << nPatches
                    << abort(FatalIOError);
            }
        }
    }

    if (nPolyPatches != nPatchBoundaryModels_)
    {
        FatalIOErrorInFunction(uspBoundariesDict_)
            << " Number of poly-patches = "  << nPolyPatches
            << " are not equal to the number of patch "
            << "models = " << nPatchBoundaryModels_
            << abort(FatalIOError);
    }
}


void uspBoundaries::setInitialConfig()
{
    for (auto& model : patchBoundaryModels_)
    {
        model->initialConfiguration();
    }

    for (auto& model : cyclicBoundaryModels_)
    {
        model->initialConfiguration();
    }

    for (auto& model : generalBoundaryModels_)
    {
        model->initialConfiguration();
    }
}


void uspBoundaries::calculateProps()
{
    for (auto& model : patchBoundaryModels_)
    {
        model->calculateProperties();
    }

    for (auto& model : cyclicBoundaryModels_)
    {
        model->calculateProperties();
    }

    for (auto& model : generalBoundaryModels_)
    {
        model->calculateProperties();
    }
}


void uspBoundaries::controlBeforeMove()
{
    // Impose model after calculation of forces
    for (auto& model: generalBoundaryModels_)
    {
        model->controlParcelsBeforeMove();
    }
}


void uspBoundaries::controlBeforeCollisions()
{
    for (auto& model : generalBoundaryModels_)
    {
        model->controlParcelsBeforeCollisions();
    }
}


void uspBoundaries::controlAfterCollisions()
{
    for (auto& model : generalBoundaryModels_)
    {
        model->controlParcelsAfterCollisions();
    }
}


void uspBoundaries::outputResults()
{
    if (!time_.writeTime())
    {
        return;
    }

    {
        List<fileName> timePathNames(pBFixedPathNames_.size());

        if (nPatchBoundaryModels_ > 0)
        {
            // directory: case/<timeDir>/uniform
            fileName uniformTimePath
                (time_.path()/time_.timeName()/"uniform");

            if (!isDir(uniformTimePath))
            {
                mkDir(uniformTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries
            fileName boundariesTimePath(uniformTimePath/"boundaries");

            if (!isDir(boundariesTimePath))
            {
                mkDir(boundariesTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries/usp
            fileName uspTimePath(boundariesTimePath/"usp");

            if (!isDir(uspTimePath))
            {
                mkDir(uspTimePath);
            }

            // directory: case/<timeDir>/uniform
                    // boundaries/usp/patchBoundaryModels
            fileName uspPatchBoundaryModelsTimePath
                (uspTimePath/"patchBoundaryModels");

            if (!isDir(uspPatchBoundaryModelsTimePath))
            {
                mkDir(uspPatchBoundaryModelsTimePath);
            }

            forAll(patchBoundaryModels_, p)
            {
                if
                (
                    patchBoundaryModels_[p]->writeInTimeDir()
                 || patchBoundaryModels_[p]->writeInCase()
                )
                {

                    // directory: case/<timeDir>/uniform/boundaries
                    // usp/patchBoundaryModels/<patchBoundaryModel>
                    fileName pBTimePath
                    (
                        uspPatchBoundaryModelsTimePath
                       /patchBoundaryNames_[p]
                    );

                    if (!isDir(pBTimePath))
                    {
                        mkDir(pBTimePath);
                    }
                    // creating directory for different zones but of the
                    // same model
                    const word& patchName =
                        patchBoundaryModels_[p]->patchName();

                    // directory: case/<timeDir>/uniform/
                                // boundaries/usp/
                                // patchBoundaryModels/
                                //<patchBoundaryModel>/<patchName>
                    fileName patchTimePath
                    (
                        pBTimePath/patchName
                    );

                    if (!isDir(patchTimePath))
                    {
                        mkDir(patchTimePath);
                    }

                    timePathNames[p] = patchTimePath;

                    patchBoundaryModels_[p]->output
                        (pBFixedPathNames_[p], timePathNames[p]);
                }
            }
        }
    }

    {
        List<fileName> timePathNames(cMFixedPathNames_.size());

        if (nCyclicBoundaryModels_ > 0)
        {
            // directory: case/<timeDir>/uniform
            fileName uniformTimePath
                (time_.path()/time_.timeName()/"uniform");

            if (!isDir(uniformTimePath))
            {
                mkDir(uniformTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries
            fileName boundariesTimePath(uniformTimePath/"boundaries");

            if (!isDir(boundariesTimePath))
            {
                mkDir(boundariesTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries/usp
            fileName uspTimePath(boundariesTimePath/"usp");

            if (!isDir(uspTimePath))
            {
                mkDir(uspTimePath);
            }

            // directory: case/<timeDir>/uniform
            // boundaries/usp/cyclicBoundaryModels
            fileName uspCyclicBoundaryModelsTimePath
                (uspTimePath/"cyclicBoundaryModels");

            if (!isDir(uspCyclicBoundaryModelsTimePath))
            {
                mkDir(uspCyclicBoundaryModelsTimePath);
            }

            forAll(cyclicBoundaryModels_, c)
            {
                if
                (
                    cyclicBoundaryModels_[c]->writeInTimeDir()
                 || cyclicBoundaryModels_[c]->writeInCase()
                )
                {
                    // directory: case/<timeDir>/uniform/boundaries
                    // usp/cyclicBoundaryModels/<cyclicBoundaryModel>
                    fileName cMTimePath
                    (
                        uspCyclicBoundaryModelsTimePath
                       /cyclicBoundaryNames_[c]
                    );

                    if (!isDir(cMTimePath))
                    {
                        mkDir(cMTimePath);
                    }

                    // creating directory for different zones but
                    // of the same model
                    const word& patchName =
                        cyclicBoundaryModels_[c]->patchName();

                    // directory: case/<timeDir>/uniform/boundaries
                    // usp/cyclicBoundaryModels
                    //<cyclicBoundaryModel>/<patchName>
                    fileName patchTimePath(cMTimePath/patchName);

                    if (!isDir(patchTimePath))
                    {
                        mkDir(patchTimePath);
                    }

                    timePathNames[c] = patchTimePath;

                    cyclicBoundaryModels_[c]->output
                    (
                        cMFixedPathNames_[c],
                        timePathNames[c]
                    );
                }
            }
        }
    }

    {
        List<fileName> timePathNames(gMFixedPathNames_.size());

        if (nGeneralBoundaryModels_ > 0)
        {
            // directory: case/<timeDir>/uniform
            fileName uniformTimePath
            (
                time_.path()/time_.timeName()/"uniform"
            );

            if (!isDir(uniformTimePath))
            {
                mkDir(uniformTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries
            fileName boundariesTimePath(uniformTimePath/"boundaries");

            if (!isDir(boundariesTimePath))
            {
                mkDir(boundariesTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries/usp
            fileName uspTimePath(boundariesTimePath/"usp");

            if (!isDir(uspTimePath))
            {
                mkDir(uspTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries
            // usp/generalBoundaryModels
            fileName uspGeneralBoundaryModelsTimePath
            (
                uspTimePath/"generalBoundaryModels"
            );

            if (!isDir(uspGeneralBoundaryModelsTimePath))
            {
                mkDir(uspGeneralBoundaryModelsTimePath);
            }

            forAll(generalBoundaryModels_, g)
            {
                if
                (
                    generalBoundaryModels_[g]->writeInTimeDir()
                 || generalBoundaryModels_[g]->writeInCase()
                )
                {
                    // directory: case/<timeDir>/uniform/boundaries
                    // usp/generalBoundaryModels/<generalBoundaryModel>
                    fileName gMTimePath
                    (
                        uspGeneralBoundaryModelsTimePath
                        /generalBoundaryNames_[g]          
                    );

                    if (!isDir(gMTimePath))
                    {
                        mkDir(gMTimePath);
                    }

                    // Creating directory for different zones but
                    // of the same model
                    const word& patchName =
                        generalBoundaryModels_[g]->patchName();

                    // directory: case/<timeDir>/uniform/boundaries
                    // usp/generalBoundaryModels
                    // <generalBoundaryModel>/<patchName>
                    fileName patchTimePath(gMTimePath/patchName);

                    if (!isDir(patchTimePath))
                    {
                        mkDir(patchTimePath);
                    }

                    timePathNames[g] = patchTimePath;

                    generalBoundaryModels_[g]->output
                    (
                        gMFixedPathNames_[g], timePathNames[g]
                    );
                }
            }
        }
    }

    {
        patchBoundaryList_.clear();

        patchBoundaryList_ =
            uspBoundariesDict_.lookup("uspPatchBoundaries");

        forAll(patchBoundaryModels_, p)
        {
            const entry& boundaryI = patchBoundaryList_[p];
            const dictionary& boundaryIDict = boundaryI.dict();

            patchBoundaryModels_[p]->updateProperties(boundaryIDict);
        }
    }

    {
        cyclicBoundaryList_.clear();

        cyclicBoundaryList_ =
            uspBoundariesDict_.lookup("uspCyclicBoundaries");

        forAll(cyclicBoundaryModels_, c)
        {
            const entry& boundaryI = cyclicBoundaryList_[c];
            const dictionary& boundaryIDict = boundaryI.dict();

            cyclicBoundaryModels_[c]->updateProperties(boundaryIDict);
        }
    }

    {
        generalBoundaryList_.clear();

        generalBoundaryList_ =
            uspBoundariesDict_.lookup("uspGeneralBoundaries");

        forAll(generalBoundaryModels_, g)
        {
            const entry& boundaryI = generalBoundaryList_[g];
            const dictionary& boundaryIDict = boundaryI.dict();

            generalBoundaryModels_[g]->updateProperties(boundaryIDict);
        }
    }
}


label uspBoundaries::nPatchBoundaryModels() const
{
    return nPatchBoundaryModels_;
}


label uspBoundaries::nCyclicBoundaryModels() const
{
    return nCyclicBoundaryModels_;
}


label uspBoundaries::nGeneralBoundaryModels() const
{
    return nGeneralBoundaryModels_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
