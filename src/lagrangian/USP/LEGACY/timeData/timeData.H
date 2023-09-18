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

Class
    timeData

Description
    Used as a time class framework for the controllers.


SourceFiles
    timeData.C

\*---------------------------------------------------------------------------*/

#ifndef timeData_H
#define timeData_H

#include "Time.H"
#include "scalarField.H"
#include "timeInterval.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                          Class timeData Declaration
\*---------------------------------------------------------------------------*/

class timeData
{
    // Private data

        const Time& time_;

        dictionary timeDict_;

        scalar writeInterval_;
        label writeIntSteps_;

        bool resetFieldsAtOutput_;
        label resetIndex_;

        timeInterval mdTime_;
        timeInterval samplingTime_;
        timeInterval averagingTime_;
        timeInterval controlTime_;

        dimensionedScalar nAvTimeSteps_;

        // this is the ratio of number of averaging steps to the number
        // of control steps (or their time-intervals)
        scalar nControlSteps_;

        label totalNSampSteps_;
        label totalNAvSteps_;
        label totalNContSteps_;

        label controlTimeIndex_;
        label averagingTimeIndex_;

        bool decoupledFromWriteInterval_;

public:

    // Constructors

        //- constructor from Time
        timeData(const Time& t);

        //- Construct from Time and timeDict
        timeData(const Time& t, const dictionary& timeDict);


    // Destructor
    ~timeData() = default;


    // Member functions

        void setInitialData();

        void checkAndModifyTimeProperties();

        void setTimeData(const dictionary& timeDict);


        //-return constant access to time
        const Time& time() const;

        //- write interval
        scalar writeInterval() const;
        label writeIntervalSteps() const;


        bool samplingTime() const;
        bool averagingTime() const;
        bool controlTime() const;

        //- return the number of ...
        label nSamples() const;
        label nControls() const;
        label nAverages() const;

        const dimensionedScalar& nAvTimeSteps() const;
        scalar nControlSteps() const;

        // for resetting
        scalar nAveragingTimeSteps();


        label totalNSampSteps() const;
        label totalNAvSteps() const;
        label totalNContSteps() const;

        label controlTimeIndex() const;
        label averagingTimeIndex() const;

        const timeInterval& mdTimeInterval() const;
        const timeInterval& sampleTimeInterval() const;
        const timeInterval& averageTimeInterval() const;
        const timeInterval& controlTimeInterval() const;

        timeInterval& controlTimeInterval();
        timeInterval& averageTimeInterval();


        bool resetFieldsAtOutput() const;
        bool& resetFieldsAtOutput();

        scalarField controlTimes();
        scalarField averagingTimes();


    // Member operators

        //- Prefix increment
        timeData& operator++();

        //- Postfix increment
        timeData& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
