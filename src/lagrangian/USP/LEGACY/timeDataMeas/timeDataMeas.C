/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "timeDataMeas.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeDataMeas::timeDataMeas
(
    const Time& t
)
:
    time_(t),
    timeDict_(),
    timeMeasOption_(),
    writeInterval_(t.controlDict().get<scalar>("writeInterval")),
    writeIntSteps_(label((writeInterval_/t.deltaT().value())  + 0.5)),
    resetFieldsAtOutput_(true),
    resetFieldsAtOutputUntilTime_(VGREAT),
    resetIndex_(0),
    mdTime_(1),
    samplingTime_(),
    averagingTime_(),
    writeTime_(),
    nAvTimeSteps_("nAvTimeSteps_", dimless, 0.0),
    totalNSampSteps_(0),
    totalNAvSteps_(0),
    averagingTimeIndex_(0),
    samplingTimeIndex_(0)
{}


//- Construct from Time and timeDict
timeDataMeas::timeDataMeas
(
    const Time& t,
    const dictionary& timeDict
)
:
    time_(t),
    timeDict_(timeDict),
    timeMeasOption_(timeDict_.get<word>("timeOption")),
    writeInterval_(t.controlDict().get<scalar>("writeInterval")),
    writeIntSteps_(label((writeInterval_/t.deltaT().value()) + 0.5)),
    resetFieldsAtOutput_(true),
    resetFieldsAtOutputUntilTime_(VGREAT),
    resetIndex_(0),
    mdTime_(1),
    samplingTime_(),
    averagingTime_(),
    writeTime_(),
    nAvTimeSteps_("nAvTimeSteps_", dimless, 0.0),
    totalNSampSteps_(0),
    totalNAvSteps_(0),
    averagingTimeIndex_(0),
    samplingTimeIndex_(0)
{
    setInitialData();
}


void timeDataMeas::setInitialData()
{
    Info << nl << "TimeData Statistics: " << endl;

    scalar deltaTMD = time_.controlDict().get<scalar>("deltaT");

    mdTime_.deltaT() = deltaTMD;

    timeDict_.readIfPresent<bool>("resetAtOutput", resetFieldsAtOutput_);
    timeDict_.readIfPresent<scalar>("resetAtOutputUntilTime", resetFieldsAtOutputUntilTime_);

    if (timeMeasOption_ == "write")
    {
        samplingTime_.nSteps() = 1;
        averagingTime_.nSteps() = writeIntSteps_;
    }
    else if (timeMeasOption_ == "decoupledWrite")
    {
        samplingTime_.nSteps() = 1;

        IOdictionary couplingDict
        (
            IOobject
            (
                "couplingDict",
                time().system(),
                time().db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dictionary timeCouplingDict
                (couplingDict.subDict("timeCouplingProperties"));

        const scalar couplingTimeInterval =
            timeCouplingDict.get<scalar>("couplingTimeInterval");
        const scalar molecularOnInterval =
            timeCouplingDict.get<scalar>("molecularOnInterval");

        scalar ratio = writeInterval_/couplingTimeInterval;
        label nAveragingSteps = label((ratio*molecularOnInterval)/deltaTMD);

        averagingTime_.nSteps() = nAveragingSteps;
    }
    else if (timeMeasOption_ == "general")
    {
        const label nSamples = timeDict_.get<label>("nSamples");
        const label nAverages = timeDict_.get<label>("nAverages");

        samplingTime_.nSteps() = nSamples;
        averagingTime_.nSteps() = nAverages;

        checkAndModifyTimeProperties();
    }
    else if (timeMeasOption_ == "coupling")
    {
        const label nSamplesDict = timeDict_.get<scalar>("nSamples");

        IOdictionary couplingDict
        (
            IOobject
            (
                "couplingDict",
                time().system(),
                time().db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dictionary timeCouplingDict
        (
            couplingDict.subDict("timeCouplingProperties")
        );

        const scalar averagingTime =
            timeCouplingDict.get<scalar>("molecularCouplingTimeInterval");

        samplingTime_.nSteps() = nSamplesDict;
        averagingTime_.nSteps() = label((averagingTime/deltaTMD) + 0.5);

        Info << "Averaging time: " << averagingTime << ", deltaTMD: "
        << deltaTMD << averagingTime_.nSteps() << endl;

        //- test properties

        bool changedProperties = false;

        const label nAverages = averagingTime_.nSteps();
        label& nSamples = samplingTime_.nSteps();

        if (nSamples < 1)
        {
            nSamples = 1;
            changedProperties = true;
        }
        else
        {
            if (nSamples > nAverages)
            {
                nSamples = nAverages;
                changedProperties = true;
            }
            else
            {
                while ((nAverages % nSamples) != 0)
                {
                    nSamples--;
                    changedProperties = true;
                }
            }
        }

        if (changedProperties)
        {
            Info << "initial nSamples: " << nSamplesDict
                << ", modified nSamples: " << samplingTime_.nSteps() << endl;

            FatalErrorInFunction
                << "Time properties are inconsistent."
                << " Check and change them appropriately from the time "
                << "dictionary"
                << nl
                << exit(FatalError);
        }
    }

    else if (timeMeasOption_ == "decoupled")
    {
        const label nSamples = timeDict_.get<label>("nSamples");

        IOdictionary couplingDict
        (
            IOobject
            (
                "couplingDict",
                time().system(),
                time().db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dictionary timeCouplingDict
        (
            couplingDict.subDict("timeCouplingProperties")
        );

        const scalar averagingTime =
            timeCouplingDict.get<scalar>("molecularOnInterval");

        samplingTime_.nSteps() = nSamples;
        averagingTime_.nSteps() = label(averagingTime/deltaTMD);

        checkAndModifyTimeProperties();
    }
    else if (timeMeasOption_ == "decoupledFromWriteInterval")
    {
        const label nSamplesDict = timeDict_.get<label>("nSamples");
        const label nAveragesDict = timeDict_.get<label>("nAverages");

        samplingTime_.nSteps() = nSamplesDict;
        averagingTime_.nSteps() = nAveragesDict;

        //- test properties

        bool changedProperties = false;

        const label nAverages = averagingTime_.nSteps();
        label& nSamples = samplingTime_.nSteps();

        if (nSamples < 1)
        {
            nSamples = 1;
            changedProperties = true;
        }
        else
        {
            if (nSamples > nAverages)
            {
                nSamples = nAverages;
                changedProperties = true;
            }
            else
            {
                while ((nAverages % nSamples) != 0)
                {
                    nSamples--;
                    changedProperties = true;
                }
            }
        }

        if (changedProperties)
        {
            Info << "initial nSamples: " << nSamplesDict
                << ", modified nSamples: " << samplingTime_.nSteps() << endl;

            FatalErrorInFunction
                << "Time properties are inconsistent."
                << " Check and change them appropriately "
                << "from the time dictionary"
                << nl
                << exit(FatalError);
        }
    }

    else if (timeMeasOption_ == "numberWrite")
    {
        const label noWriteIntervals =
            timeDict_.get<label>("noOfWriteIntervals");
        const label nSamplesDict = timeDict_.get<label>("nSamples");

        averagingTime_.nSteps() =
            (scalar(noWriteIntervals)*writeInterval_)/mdTime_.deltaT();
        samplingTime_.nSteps() = nSamplesDict;
        writeTime_.nSteps() = averagingTime_.nSteps();

        //- test properties

        bool changedProperties = false;

        const label nAverages = averagingTime_.nSteps();
        label& nSamples = samplingTime_.nSteps();

        if (nSamples < 1)
        {
            nSamples = 1;
            changedProperties = true;
        }
        else
        {
            if (nSamples > nAverages)
            {
                nSamples = nAverages;
                changedProperties = true;
            }
            else
            {
                while ((nAverages % nSamples) != 0)
                {
                    nSamples--;
                    changedProperties = true;
                }
            }
        }

        if (changedProperties)
        {
            Info << "initial nSamples: " << nSamplesDict
                << ", modified nSamples: " << samplingTime_.nSteps() << endl;

            FatalErrorInFunction
                << "Time properties are inconsistent."
                << " Check and change them appropriately from the time "
                << "dictionary"
                << nl
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "timeOption: \"" << timeMeasOption_
            << "\" is not one of the available options." << nl
            << exit(FatalError);
    }

    samplingTime_.deltaT() = deltaTMD * scalar(samplingTime_.nSteps());
    averagingTime_.deltaT() = deltaTMD * scalar(averagingTime_.nSteps());

    Info << " measurement option: " << timeMeasOption_ << endl;
    Info << " nSamples: " << samplingTime_.nSteps()
         << ", time interval: " << samplingTime_.deltaT()
         << endl;

    Info << " nAverages: " << averagingTime_.nSteps()
         << ", time interval: " << averagingTime_.deltaT()
         << endl;


    const scalar endTime = time_.endTime().value();
    const scalar startTime = time_.startTime().value();

    totalNAvSteps_ = label (((endTime-startTime) / averagingTime_.deltaT())
                                                                        + 0.5);

    totalNSampSteps_ = label(((endTime-startTime) / samplingTime_.deltaT())
                                                                        + 0.5);

    Info << " total no. of sampling steps: " << totalNSampSteps_ << endl;
    Info << " total no. of averaging Steps: " << totalNAvSteps_ << endl;
    Info << nl << endl;

    nAvTimeSteps_.value() = scalar(averagingTime_.nSteps());
}


void timeDataMeas::checkAndModifyTimeProperties()
{
    //- checking

    bool changedProperties = false;

    const label nAveragesOriginal = averagingTime_.nSteps();
    const label nSamplesOriginal = samplingTime_.nSteps();

    // - 1. averaging time

    label& nAverages = averagingTime_.nSteps();

    if (nAverages < 1)
    {
        nAverages = 1;
        changedProperties = true;
    }
    else
    {
        if (nAverages > writeIntSteps_)
        {
            nAverages = writeIntSteps_;
            changedProperties = true;
        }
        else
        {
            while ((writeIntSteps_ % nAverages) != 0)
            {
                nAverages++;
                changedProperties = true;
            }
        }
    }

    // - 2. sampling time
    label& nSamples = samplingTime_.nSteps();

    if (nSamples < 1)
    {
        nSamples = 1;
        changedProperties = true;
    }
    else
    {
        if (nSamples > nAverages)
        {
            nSamples = nAverages;
            changedProperties = true;
        }
        else
        {
            while ((nAverages % nSamples) != 0)
            {
                nSamples--;
                changedProperties = true;
            }
        }
    }

    if (changedProperties)
    {
        Info << "initial nSamples: " << nSamplesOriginal
             << ", modified nSamples: " << samplingTime_.nSteps() << endl;

        Info << "initial nAverages: " << nAveragesOriginal
             << ", modified nAverages: " << averagingTime_.nSteps() << endl;

        FatalErrorInFunction
            << "Time properties are inconsistent."
            << " Check and change them appropriately from the time dictionary"
            << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

timeDataMeas::~timeDataMeas()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeDataMeas::setTimeData(const dictionary& timeDict)
{
    const label nSamples(timeDict.get<label>("nSamples"));
    const label nAverages(timeDict.get<label>("nAverages"));

    samplingTime_.nSteps() = nSamples;
    averagingTime_.nSteps() = nAverages;

    setInitialData();
}


const Time& timeDataMeas::time() const
{
    return time_;
}


scalar timeDataMeas::writeInterval() const
{
    return writeInterval_;
}


label timeDataMeas::writeIntervalSteps() const
{
    return writeIntSteps_;
}


bool timeDataMeas::samplingTime() const
{
    return samplingTime_.endTime();
}


bool timeDataMeas::averagingTime() const
{
    return averagingTime_.endTime();
}


bool timeDataMeas::writeTime() const
{
    return writeTime_.endTime();
}


label timeDataMeas::nSamples() const
{
    return samplingTime_.nSteps();
}


label timeDataMeas::nAverages() const
{
    return averagingTime_.nSteps();
}


const dimensionedScalar& timeDataMeas::nAvTimeSteps() const
{
    return nAvTimeSteps_;
}


// for resetting
scalar timeDataMeas::nAveragingTimeSteps()
{
    return scalar(nAvTimeSteps().value()*resetIndex_);
}


label timeDataMeas::totalNSampSteps() const
{
    return totalNSampSteps_;
}


label timeDataMeas::totalNAvSteps() const
{
    return totalNAvSteps_;
}


label timeDataMeas::averagingTimeIndex() const
{
    return averagingTimeIndex_;
}


label timeDataMeas::samplingTimeIndex() const
{
    return samplingTimeIndex_;
}


bool timeDataMeas::resetFieldsAtOutput() const
{
    return resetFieldsAtOutput_;
}


bool& timeDataMeas::resetFieldsAtOutput()
{
    return resetFieldsAtOutput_;
}

scalar timeDataMeas::resetFieldsAtOutputUntilTime() const
{
    return resetFieldsAtOutputUntilTime_;
}


scalar& timeDataMeas::resetFieldsAtOutputUntilTime()
{
    return resetFieldsAtOutputUntilTime_;
}


scalarField timeDataMeas::averagingTimes()
{
    const scalar startTime = time_.startTime().value();

    scalarField averagingTimes(totalNAvSteps_ + 1, 0.0);

    forAll(averagingTimes, tT)
    {
        averagingTimes[tT] = startTime + tT*averagingTime_.deltaT();
    }

    return averagingTimes;
}


scalarField timeDataMeas::samplingTimes()
{
    const scalar startTime = time_.startTime().value();

    scalarField samplingTimes(totalNSampSteps_ + 1, 0.0);

    forAll(samplingTimes, tT)
    {
        samplingTimes[tT] = startTime + tT*samplingTime_.deltaT();
    }

    return samplingTimes;
}


scalarField timeDataMeas::writeTimes()
{
    const scalar endTime = time_.endTime().value();
    const scalar startTime = time_.startTime().value();

    label nWriteSteps = label((endTime-startTime)/writeInterval_);

    scalarField writingTimes(nWriteSteps + 1, 0.0);

    forAll(writingTimes, tT)
    {
        writingTimes[tT] = startTime + tT*writeInterval_;
    }

    return writingTimes;
}


label timeDataMeas::writeSteps()
{
    const scalar endTime = time_.endTime().value();
    const scalar startTime = time_.startTime().value();

    label nWriteSteps = label((endTime-startTime)/writeInterval_);

    return nWriteSteps;
}


const timeInterval& timeDataMeas::mdTimeInterval() const
{
    return mdTime_;
}


const timeInterval& timeDataMeas::sampleTimeInterval() const
{
    return samplingTime_;
}


const timeInterval& timeDataMeas::averageTimeInterval() const
{
    return averagingTime_;
}


const timeInterval& timeDataMeas::writeTimeInterval() const
{
    return writeTime_;
}


const word& timeDataMeas::timeOption() const
{
    return timeMeasOption_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//- Prefix increment
timeDataMeas& timeDataMeas::operator++()
{
    samplingTime_++;
    averagingTime_++;
    writeTime_++;

//     calcPropTime_++;
    if (averagingTime_.endTime())
    {
        if (resetFieldsAtOutput_)
        {
            resetIndex_ = -1;
        }

        averagingTimeIndex_++;
        resetIndex_++;
    }
    if (samplingTime_.endTime())
    {
        samplingTimeIndex_++;
    }

    return *this;
}


//- Postfix increment
timeDataMeas& timeDataMeas::operator++(int)
{
    return operator++();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
