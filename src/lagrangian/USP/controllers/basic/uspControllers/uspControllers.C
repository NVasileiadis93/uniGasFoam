/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "uspControllers.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::uspControllers::uspControllers(const Time& t, const polyMesh& mesh)
:
    time_(t),
    uspControllersDict_
    (
        IOobject
        (
            "controllersDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    nStateControllers_(0),
    nFluxControllers_(0),
    stateControllersList_(),
    sCNames_(),
    sCIds_(),
    sCFixedPathNames_(),
    stateControllers_(),
    fluxControllersList_(),
    fCNames_(),
    fCIds_(),
    fCFixedPathNames_(),
    fluxControllers_()
{}


Foam::uspControllers::uspControllers
(
    const Time& t,
    const polyMesh& mesh,
    uspCloud& cloud
)
:
    time_(t),
    uspControllersDict_
    (
        IOobject
        (
            "controllersDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    nStateControllers_(0),
    nFluxControllers_(0),
    stateControllersList_(uspControllersDict_.lookup("uspStateControllers")),
    sCNames_(stateControllersList_.size()),
    sCIds_(stateControllersList_.size()),
    sCFixedPathNames_(stateControllersList_.size()),
    stateControllers_(stateControllersList_.size()),
    fluxControllersList_(uspControllersDict_.lookup("uspFluxControllers")),
    fCNames_(fluxControllersList_.size()),
    fCIds_(fluxControllersList_.size()),
    fCFixedPathNames_(fluxControllersList_.size()),
    fluxControllers_(fluxControllersList_.size())
{
    Info << nl << "Creating uspControllers" << nl << endl;

    // state uspControllers
    nStateControllers_ = stateControllers_.size();

    forAll(stateControllers_, sC)
    {
        const entry& uspControllersI = stateControllersList_[sC];
        const dictionary& uspControllersIDict = uspControllersI.dict();

        stateControllers_[sC] =
            uspStateController::New(cloud.mesh(), cloud, uspControllersIDict);

        sCNames_[sC] = stateControllers_[sC]->type();
        sCIds_[sC] = sC;
    }


    // flux uspControllers

    nFluxControllers_ = fluxControllers_.size();
    forAll(fluxControllers_, fC)
    {
        const entry& uspControllersI = fluxControllersList_[fC];

        const dictionary& uspControllersIDict = uspControllersI.dict();

        fluxControllers_[fC] =
            uspFluxController::New(cloud.mesh(), cloud, uspControllersIDict);

        fCNames_[fC] = fluxControllers_[fC]->type();
        fCIds_[fC] = fC;
    }

    // creating directories for state controllers
    if (nStateControllers_ > 0)
    {
        // directory: case/controllers
        fileName controllersPath(time_.path()/"controllers");

        if (!isDir(controllersPath))
        {
            mkDir(controllersPath);
        }

        // directory: case/controllers/usp
        fileName uspControllersPath(controllersPath/"usp");

        if (!isDir(uspControllersPath))
        {
            mkDir(uspControllersPath);
        }

        // directory: case/controllers/usp/stateControllers
        fileName stateControllersPath(uspControllersPath/"stateControllers");

        if (!isDir(stateControllersPath))
        {
            mkDir(stateControllersPath);
        }

        forAll(stateControllers_, sC)
        {
            if (stateControllers_[sC]->writeInCase())
            {
                // directory:
                // case/controllers/usp/stateControllers/<stateControllerModel>
                fileName stateControllerPath(stateControllersPath/sCNames_[sC]);

                if (!isDir(stateControllerPath))
                {
                    mkDir(stateControllerPath);
                }

                const word& regionName = stateControllers_[sC]->regionName();

                // directory:
                // case/controllers/usp/stateControllers/
                // <stateControllerModel>/<cellZoneName>
                fileName zonePath(stateControllerPath/regionName);

                if (!isDir(zonePath))
                {
                    mkDir(zonePath);
                }

                sCFixedPathNames_[sC] = zonePath;
            }
        }
    }

    // creating directories for flux controllers
    if (nFluxControllers_ > 0)
    {
        // directory: case/controllers
        fileName controllersPath(time_.path()/"controllers");

        if (!isDir(controllersPath))
        {
            mkDir(controllersPath);
        }

        // directory: case/controllers/usp
        fileName uspControllersPath(time_.path()/"usp");

        if (!isDir(uspControllersPath))
        {
            mkDir(uspControllersPath);
        }

        // directory: case/controllers/usp/fluxControllers
        fileName fluxControllersPath(uspControllersPath/"fluxControllers");

        if (!isDir(fluxControllersPath))
        {
            mkDir(fluxControllersPath);
        }

        forAll(fluxControllers_, fC)
        {
            if (fluxControllers_[fC]->writeInCase())
            {
                // directory:
                // case/controllers/usp/fluxControllers/<fluxControllerModel>
                fileName fluxControllerPath(fluxControllersPath/fCNames_[fC]);

                if (!isDir(fluxControllerPath))
                {
                    mkDir(fluxControllerPath);
                }

                const word& regionName = fluxControllers_[fC]->regionName();

                // directory:
                // case/controllers/usp/fluxControllers/
                // <fluxControllerModel>/<faceZoneName>
                fileName zonePath(fluxControllerPath/regionName);

                if (!isDir(zonePath))
                {
                    mkDir(zonePath);
                }

                fCFixedPathNames_[fC] = zonePath;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspControllers::initialConfig()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->initialConfiguration();
    }

    forAll(fluxControllers_, fC)
    {
        fluxControllers_[fC]->initialConfiguration();
    }
}


void Foam::uspControllers::controlBeforeMove()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlParcelsBeforeMove();
    }
}


void Foam::uspControllers::controlBeforeCollisions()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlParcelsBeforeCollisions();
    }
}


void Foam::uspControllers::controlAfterCollisions()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlParcelsAfterCollisions();
    }
}


void Foam::uspControllers::calculateProps()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->calculateProperties();
    }

    forAll(fluxControllers_, fC)
    {
        fluxControllers_[fC]->calculateProperties();
    }
}


void Foam::uspControllers::updateTimeInfo()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->updateTime();
    }

    forAll(fluxControllers_, fC)
    {
        fluxControllers_[fC]->updateTime();
    }
}


void Foam::uspControllers::outputResults()
{
    if (!time_.writeTime())
    {
        return;
    }

    // Creating a set of directories in the current time directory
    {
        List<fileName> timePathNames(sCFixedPathNames_.size());

        if (nStateControllers_ > 0)
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

            // directory: case/<timeDir>/uniform/controllers
            fileName controllersTimePath(uniformTimePath/"controllers");

            if (!isDir(controllersTimePath))
            {
                mkDir(controllersTimePath);
            }

            // directory: case/<timeDir>/uniform/controllers/usp
            fileName uspTimePath(controllersTimePath/"usp");

            if (!isDir(uspTimePath))
            {
                mkDir(uspTimePath);
            }

            // directory: case/<timeDir>/uniform/controllers/usp/
            fileName uspStateControllersTimePath
            (
                uspTimePath/"stateControllers"
            );

            if (!isDir(uspStateControllersTimePath))
            {
                mkDir(uspStateControllersTimePath);
            }

            forAll(stateControllers_, sC)
            {
                if (stateControllers_[sC]->writeInTimeDir())
                {
                    // directory:
                    // case/<timeDir>/uniform/controllers/
                    // usp/<stateControllerModel>
                    fileName sCTimePath =
                        uspStateControllersTimePath/sCNames_[sC];

                    if (!isDir(sCTimePath))
                    {
                        mkDir(sCTimePath);
                    }

                    // creating directory for different zones but
                    // of the same model
                    const word& regionName =
                        stateControllers_[sC]->regionName();

                    // directory: case/<timeDir>/uniform/controllers/
                    // usp/<stateControllerModel>/
                    // <cellZoneName>
                    fileName zoneTimePath(sCTimePath/regionName);

                    if (!isDir(zoneTimePath))
                    {
                        mkDir(zoneTimePath);
                    }

                    timePathNames[sC] = zoneTimePath;
                }
            }
        }

        // Write data (do not comment this out)
        forAll(stateControllers_, sC)
        {
            stateControllers_[sC]->output
            (
                sCFixedPathNames_[sC],
                timePathNames[sC]
            );
        }
    }

    {
        List<fileName> timePathNames(fCFixedPathNames_.size());

        if (nFluxControllers_ > 0)
        {
            // directory: case/<timeDir>/uniform
            fileName uniformTimePath =
                time_.path()/time_.timeName()/"uniform";

            if (!isDir(uniformTimePath))
            {
                mkDir(uniformTimePath);
            }

            // directory: case/<timeDir>/uniform/controllers
            fileName controllersTimePath(uniformTimePath/"controllers");

            if (!isDir(controllersTimePath))
            {
                mkDir(controllersTimePath);
            }

            // directory: case/<timeDir>/uniform/controllers/usp
            fileName uspTimePath(controllersTimePath/"usp");

            if (!isDir(uspTimePath))
            {
                mkDir(uspTimePath);
            }

            // directory: case/<timeDir>/uniform/fluxControllers
            fileName uspControllersTimePath
            (
                uspTimePath/"fluxControllers"
            );

            if (!isDir(uspControllersTimePath))
            {
                mkDir(uspControllersTimePath);
            }

            forAll(fluxControllers_, fC)
            {
                if (fluxControllers_[fC]->writeInTimeDir())
                {
                    // directory: case/<timeDir>/uniform/controllers/
                    // usp/<fluxControllerModel>
                    fileName fCTimePath
                    (
                        uspControllersTimePath/fCNames_[fC]
                    );

                    if (!isDir(fCTimePath))
                    {
                        mkDir(fCTimePath);
                    }

                    const word& regionName =
                        fluxControllers_[fC]->regionName();

                    // directory: case/<timeDir>/uniform/controllers/
                    // usp/<fluxControllerModel>
                    // <faceZoneName>
                    fileName zoneTimePath(fCTimePath/regionName);

                    if (!isDir(zoneTimePath))
                    {
                        mkDir(zoneTimePath);
                    }

                    timePathNames[fC] = zoneTimePath;
                }
            }
        }

        // Write data (do not comment this out)
        forAll(fluxControllers_, fC)
        {
            fluxControllers_[fC]->output
            (
                fCFixedPathNames_[fC],
                timePathNames[fC]
            );
        }
    }

    // Re-read dictionaries for modified properties (run-time selection)
    {
        stateControllersList_.clear();

        stateControllersList_ =
            uspControllersDict_.lookup("uspStateControllers");

        forAll(stateControllers_, sC)
        {
            const entry& uspControllersI = stateControllersList_[sC];
            const dictionary& uspControllersIDict =
                uspControllersI.dict();

            stateControllers_[sC]->updateProperties(uspControllersIDict);
        }
    }

    {
        fluxControllersList_.clear();

        fluxControllersList_ =
            uspControllersDict_.lookup("uspFluxControllers");

        forAll(fluxControllers_, fC)
        {
            const entry& uspControllersI = fluxControllersList_[fC];
            const dictionary& uspControllersIDict =
                uspControllersI.dict();

            fluxControllers_[fC]->updateProperties(uspControllersIDict);
        }
    }
}


Foam::label Foam::uspControllers::nStateControllers() const
{
    return nStateControllers_;
}


// ************************************************************************* //
