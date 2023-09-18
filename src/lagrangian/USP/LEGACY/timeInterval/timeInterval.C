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

#include "timeInterval.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null Constructor
timeInterval::timeInterval()
:
    nSteps_(-1),
    deltaT_(0.0),
    timeIndex_(0),
    endTime_(false)
{}


// Construct from components
timeInterval::timeInterval
(
    const label nSteps
)
:
    nSteps_(nSteps),
    deltaT_(0.0),
    timeIndex_(0),
    endTime_(false)
{}


timeInterval::timeInterval
(
    const label nSteps,
    const scalar deltaT
)
:
    nSteps_(nSteps),
    deltaT_(deltaT),
    timeIndex_(0),
    endTime_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

timeInterval::~timeInterval()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool timeInterval::endTime() const
{
    return endTime_;
}


bool& timeInterval::endTime()
{
    return endTime_;
}


label timeInterval::nSteps() const
{
    return nSteps_;
}


label& timeInterval::nSteps()
{
    return nSteps_;
}


label timeInterval::timeIndex() const
{
    return timeIndex_;
}


label& timeInterval::timeIndex()
{
    return timeIndex_;
}


scalar timeInterval::deltaT() const
{
    return deltaT_;
}


scalar& timeInterval::deltaT()
{
    return deltaT_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//- Prefix increment
timeInterval& timeInterval::operator++()
{
    endTime_ = false;
    timeIndex_++;

    if (!(timeIndex_ < nSteps_))
    {
        endTime_ = true;
        timeIndex_ = 0;
    }

    return *this;
}


//- Postfix increment
timeInterval& timeInterval::operator++(int)
{
    return operator++();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
