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

#include "uspStateController.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspStateController, 0);
defineRunTimeSelectionTable(uspStateController, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspStateController::uspStateController
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspControllerBase(mesh, cloud, dict),
    timePeriod_(timeDict_.get<scalar>("initialTimePeriod")),
    initialTime_(mesh_.time().startTime().value())
{
    const cellZoneMesh& cellZones = mesh_.cellZones();
    regionId_ = cellZones.findZoneID(regionName_);

    if (regionId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find cellZone: " << regionName_ << nl << " in: "
            << mesh_.time().system()/"controllersDict"
            << exit(FatalError);
    }

    const scalar avTimeInterval = timeData_.averageTimeInterval().deltaT();

    if ((timePeriod_ < avTimeInterval) && (timePeriod_ > 0.0))
    {
        timePeriod_ = avTimeInterval;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::uspStateController> Foam::uspStateController::New
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("stateControllerModel"));

    Info<< "Selecting uspStateController " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "uspStateController",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<uspStateController>(cstrIter()(mesh, cloud, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspStateController::updateTime()
{
    uspControllerBase::updateTime();

    const scalar t = mesh_.time().timeOutputValue();

    if ((t - initialTime_) < timePeriod_)
    {
        timeData_.controlTimeInterval().endTime() = false;
    }
}


void Foam::uspStateController::updateProperties(const dictionary& dict)
{
    uspControllerBase::updateProperties(dict);

    timeDict_ = controllerDict_.subDict("timeProperties");
    timeDict_.readIfPresent("resetAtOutput", timeData_.resetFieldsAtOutput());
}


const Foam::labelList& Foam::uspStateController::controlZone() const
{
    return mesh_.cellZones()[regionId_];
}


// ************************************************************************* //
