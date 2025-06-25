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

#include "uniGasBoundaries.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::uniGasBoundaries::dictName("boundariesDict");


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

uniGasBoundaries::uniGasBoundaries(const Time& t, const polyMesh& mesh)
:
    time_(t),
    uniGasBoundariesDict_
    (
        IOobject
        (
            uniGasBoundaries::dictName,
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


uniGasBoundaries::uniGasBoundaries
(
    const Time& t,
    const polyMesh& mesh,
    uniGasCloud& cloud
)
:
    time_(t),
    uniGasBoundariesDict_
    (
        IOobject
        (
            uniGasBoundaries::dictName,
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    nPatchBoundaryModels_(0),
    nCyclicBoundaryModels_(0),
    nGeneralBoundaryModels_(0),

    patchBoundaryList_(uniGasBoundariesDict_.lookup("uniGasPatchBoundaries")),
    patchBoundaryNames_(patchBoundaryList_.size()),
    patchBoundaryIds_(patchBoundaryList_.size()),
    pBFixedPathNames_(patchBoundaryList_.size()),
    patchBoundaryModels_(patchBoundaryList_.size()),
    patchToModelId_(mesh.boundaryMesh().size(), -1),

    cyclicBoundaryList_(uniGasBoundariesDict_.lookup("uniGasCyclicBoundaries")),
    cyclicBoundaryNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryIds_(cyclicBoundaryList_.size()),
    cMFixedPathNames_(cyclicBoundaryList_.size()),
    cyclicBoundaryModels_(cyclicBoundaryList_.size()),
    cyclicBoundaryToModelId_(mesh.boundaryMesh().size(), -1),

    generalBoundaryList_(uniGasBoundariesDict_.lookup("uniGasGeneralBoundaries")),
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

            patchBoundaryModels_[p] = autoPtr<uniGasPatchBoundary>
            (
                uniGasPatchBoundary::New(mesh, cloud, boundaryIDict)
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

            cyclicBoundaryModels_[c] = autoPtr<uniGasCyclicBoundary>
            (
                uniGasCyclicBoundary::New(mesh, cloud, boundaryIDict)
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

            generalBoundaryModels_[g] = autoPtr<uniGasGeneralBoundary>
            (
                uniGasGeneralBoundary::New(mesh, cloud, boundaryIDict)
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

        // directory: case/boundaries/uniGas
        fileName uniGasBoundariesPath(boundariesPath/"uniGas");

        if (!isDir(uniGasBoundariesPath))
        {
            mkDir(uniGasBoundariesPath);
        }

        // directory: case/boundaries/uniGas/patchBoundaryModels
        fileName patchBoundaryModelsPath
        (
            uniGasBoundariesPath/"patchBoundaryModels"
        );

        if (!isDir(patchBoundaryModelsPath))
        {
            mkDir(patchBoundaryModelsPath);
        }

        forAll(patchBoundaryModels_, p)
        {
            if (patchBoundaryModels_[p]->writeInCase())
            {
                // directory: case/boundaries/uniGas/
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

                // directory: case/boundaries/uniGas/
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

        // directory: case/boundaries/uniGas
        fileName uniGasBoundariesPath(boundariesPath/"uniGas");

        if (!isDir(uniGasBoundariesPath))
        {
            mkDir(uniGasBoundariesPath);
        }

        // directory: case/boundaries/uniGas/cyclicBoundaryModels
        fileName cyclicBoundaryModelsPath
                    (uniGasBoundariesPath/"cyclicBoundaryModels");

        if (!isDir(cyclicBoundaryModelsPath))
        {
            mkDir(cyclicBoundaryModelsPath);
        }

        forAll(cyclicBoundaryModels_, c)
        {
            if (cyclicBoundaryModels_[c]->writeInCase())
            {
                // directory: case/boundaries/uniGas/
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

                // directory: case/boundaries/uniGas
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

        // directory: case/boundaries/uniGas
        fileName uniGasBoundariesPath(boundariesPath/"uniGas");

        if (!isDir(uniGasBoundariesPath))
        {
            mkDir(uniGasBoundariesPath);
        }

        // directory: case/boundaries/uniGas/cyclicBoundaryModels
        fileName generalBoundaryModelsPath
        (
            uniGasBoundariesPath/"generalBoundaryModels"
        );

        if (!isDir(generalBoundaryModelsPath))
        {
            mkDir(generalBoundaryModelsPath);
        }

        forAll(generalBoundaryModels_, g)
        {
            if (generalBoundaryModels_[g]->writeInCase())
            {
                // directory: case/boundaries/uniGas
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

                // directory: case/boundaries/uniGas
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

void uniGasBoundaries::checkCyclicBoundaryModels(const polyMesh& mesh)
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


void uniGasBoundaries::checkPatchBoundaryModels(const polyMesh& mesh)
{

    label nPolyPatches = 0;

    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        // Check that all poly-patches types are spelled correctly
        if (patch.type() == "genericPatch")
        {
                FatalIOErrorInFunction(uniGasBoundariesDict_)
                    << " Wrong type for patch: "
                    << patch.name()
                    << abort(FatalIOError);
        }

        // Check all wedge poly-patches are switched to symmetryPlane poly-patches
        if (isA<wedgePolyPatch>(patch))
        {
                FatalIOErrorInFunction(uniGasBoundariesDict_)
                    << " Wedge patch type not supported, [name: "
                    << patch.name()
                    << "]. For axisymmetric case change patch type to symmetryPlane."
                    << abort(FatalIOError);
        }

        // Check that all poly-patches defined within blockMeshDict,
        // each have one model.
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
                FatalIOErrorInFunction(uniGasBoundariesDict_)
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
        FatalIOErrorInFunction(uniGasBoundariesDict_)
            << " Number of poly-patches = "  << nPolyPatches
            << " are not equal to the number of patch "
            << "models = " << nPatchBoundaryModels_
            << abort(FatalIOError);
    }
    
}


void uniGasBoundaries::setInitialConfig()
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


void uniGasBoundaries::calculateProps()
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


void uniGasBoundaries::controlBeforeMove()
{
    // Impose model after calculation of forces
    for (auto& model: generalBoundaryModels_)
    {
        model->controlParcelsBeforeMove();
    }
}


void uniGasBoundaries::controlBeforeCollisions()
{
    for (auto& model : generalBoundaryModels_)
    {
        model->controlParcelsBeforeCollisions();
    }
}


void uniGasBoundaries::controlAfterCollisions()
{
    for (auto& model : generalBoundaryModels_)
    {
        model->controlParcelsAfterCollisions();
    }
}


void uniGasBoundaries::outputResults()
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

            // directory: case/<timeDir>/uniform/boundaries/uniGas
            fileName uniGasTimePath(boundariesTimePath/"uniGas");

            if (!isDir(uniGasTimePath))
            {
                mkDir(uniGasTimePath);
            }

            // directory: case/<timeDir>/uniform
                    // boundaries/uniGas/patchBoundaryModels
            fileName uniGasPatchBoundaryModelsTimePath
                (uniGasTimePath/"patchBoundaryModels");

            if (!isDir(uniGasPatchBoundaryModelsTimePath))
            {
                mkDir(uniGasPatchBoundaryModelsTimePath);
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
                    // uniGas/patchBoundaryModels/<patchBoundaryModel>
                    fileName pBTimePath
                    (
                        uniGasPatchBoundaryModelsTimePath
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
                                // boundaries/uniGas/
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

            // directory: case/<timeDir>/uniform/boundaries/uniGas
            fileName uniGasTimePath(boundariesTimePath/"uniGas");

            if (!isDir(uniGasTimePath))
            {
                mkDir(uniGasTimePath);
            }

            // directory: case/<timeDir>/uniform
            // boundaries/uniGas/cyclicBoundaryModels
            fileName uniGasCyclicBoundaryModelsTimePath
                (uniGasTimePath/"cyclicBoundaryModels");

            if (!isDir(uniGasCyclicBoundaryModelsTimePath))
            {
                mkDir(uniGasCyclicBoundaryModelsTimePath);
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
                    // uniGas/cyclicBoundaryModels/<cyclicBoundaryModel>
                    fileName cMTimePath
                    (
                        uniGasCyclicBoundaryModelsTimePath
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
                    // uniGas/cyclicBoundaryModels
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

            // directory: case/<timeDir>/uniform/boundaries/uniGas
            fileName uniGasTimePath(boundariesTimePath/"uniGas");

            if (!isDir(uniGasTimePath))
            {
                mkDir(uniGasTimePath);
            }

            // directory: case/<timeDir>/uniform/boundaries
            // uniGas/generalBoundaryModels
            fileName uniGasGeneralBoundaryModelsTimePath
            (
                uniGasTimePath/"generalBoundaryModels"
            );

            if (!isDir(uniGasGeneralBoundaryModelsTimePath))
            {
                mkDir(uniGasGeneralBoundaryModelsTimePath);
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
                    // uniGas/generalBoundaryModels/<generalBoundaryModel>
                    fileName gMTimePath
                    (
                        uniGasGeneralBoundaryModelsTimePath
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
                    // uniGas/generalBoundaryModels
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
            uniGasBoundariesDict_.lookup("uniGasPatchBoundaries");

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
            uniGasBoundariesDict_.lookup("uniGasCyclicBoundaries");

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
            uniGasBoundariesDict_.lookup("uniGasGeneralBoundaries");

        forAll(generalBoundaryModels_, g)
        {
            const entry& boundaryI = generalBoundaryList_[g];
            const dictionary& boundaryIDict = boundaryI.dict();

            generalBoundaryModels_[g]->updateProperties(boundaryIDict);
        }
    }
}


label uniGasBoundaries::nPatchBoundaryModels() const
{
    return nPatchBoundaryModels_;
}


label uniGasBoundaries::nCyclicBoundaryModels() const
{
    return nCyclicBoundaryModels_;
}


label uniGasBoundaries::nGeneralBoundaryModels() const
{
    return nGeneralBoundaryModels_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
