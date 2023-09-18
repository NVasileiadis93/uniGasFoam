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
    timeDataMeas

Description

SourceFiles
    timeDataMeas.C

\*---------------------------------------------------------------------------*/

#ifndef timeDataMeas_H
#define timeDataMeas_H

#include "Time.H"
#include "scalarField.H"
#include "timeInterval.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class timeDataMeas Declaration
\*---------------------------------------------------------------------------*/

class timeDataMeas
{
    // Private data

        const Time& time_;

        dictionary timeDict_;

        word timeMeasOption_;

        scalar writeInterval_;
        label writeIntSteps_;

        bool resetFieldsAtOutput_;

        scalar resetFieldsAtOutputUntilTime_;

        label resetIndex_;

        timeInterval mdTime_;
        timeInterval samplingTime_;
        timeInterval averagingTime_;
        timeInterval writeTime_;

        dimensionedScalar nAvTimeSteps_;

        label totalNSampSteps_;
        label totalNAvSteps_;

        label averagingTimeIndex_;
        label samplingTimeIndex_;
public:

    // Constructors

        //- constructor from Time
        timeDataMeas
        (
            const Time& t
        );

        //- Construct from Time and timeDict
        timeDataMeas
        (
            const Time& t,
            const dictionary& timeDict
        );


    // Destructor

        ~timeDataMeas();


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
        bool writeTime() const;

        //- return the number of ...
        label nSamples() const;
        label nAverages() const;

        const dimensionedScalar& nAvTimeSteps() const;

        bool resetFieldsAtOutput() const;
        bool& resetFieldsAtOutput();

        scalar resetFieldsAtOutputUntilTime() const;
        scalar& resetFieldsAtOutputUntilTime();

        // for resetting
        scalar nAveragingTimeSteps();



        label totalNSampSteps() const;
        label totalNAvSteps() const;

        label averagingTimeIndex() const;
        label samplingTimeIndex() const;

        const timeInterval& mdTimeInterval() const;
        const timeInterval& sampleTimeInterval() const;
        const timeInterval& averageTimeInterval() const;
        const timeInterval& writeTimeInterval() const;


        const word& timeOption() const;

        scalarField averagingTimes();
        scalarField samplingTimes();
        scalarField writeTimes();
        label writeSteps();



    // Member operators

        //- Prefix increment
        timeDataMeas& operator++();

        //- Postfix increment
        timeDataMeas& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
