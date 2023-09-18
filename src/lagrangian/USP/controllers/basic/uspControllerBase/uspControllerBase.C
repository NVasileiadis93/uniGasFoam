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

#include "uspControllerBase.H"
#include "uspCloud.H"
#include "graph.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspControllerBase::uspControllerBase
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    mesh_(mesh),
    cloud_(cloud),
    rndGen_(cloud.rndGen()),
    controllerDict_(dict.subDict("controllerProperties")),
    timeDict_(controllerDict_.subDict("timeProperties")),
    timeData_(mesh.time(), timeDict_),
    regionName_(controllerDict_.get<word>("zone")),
    regionId_(-1),
    control_(controllerDict_.get<bool>("controlSwitch")),
    readStateFromFile_(controllerDict_.get<bool>("readStateFromFile")),
    singleValueController_(false),
    density_(0.0),
    velocity_(Zero),
    temperature_(0.0),
    pressure_(0.0),
    strainRate_(Zero),
    tempGradient_(Zero),
    fieldController_(false),
    densities_(),
    velocities_(),
    temperatures_(),
    pressures_(),
    writeInTimeDir_(true),
    writeInCase_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspControllerBase::writeTimeData
(
    const fileName& pathName,
    const word& fileName,
    const List<Pair<scalar>>& data
)
{
    OFstream timeFile(pathName/fileName+".raw");

    if (timeFile.good())
    {
        timeFile << data << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Cannot open file " << timeFile.name()
            << abort(FatalError);
    }
}


void Foam::uspControllerBase::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const scalarField& yData
)
{
    fileName writeFile(pathName/nameFile);

    graph outputGraph("title", "x", "y", xData, yData);

    outputGraph.write(writeFile, "raw");
}


void Foam::uspControllerBase::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const vectorField& yData
)
{
    OFstream file(pathName/nameFile + ".xyz");

    if (file.good())
    {
        forAll(yData, n)
        {
            file<< xData[n] << "\t"
                << yData[n].x() << "\t" << yData[n].y()
                << "\t" << yData[n].z()
                << endl;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}


void Foam::uspControllerBase::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const tensorField& yData
)
{
    OFstream file(pathName/nameFile + ".xyz");

    if (file.good())
    {
        forAll(yData, n)
        {
            file<< xData[n] << "\t"
                << yData[n].xx() << "\t" << yData[n].xy()
                << "\t" << yData[n].xz() << "\t"
                << yData[n].yx() << "\t" << yData[n].yy()
                << "\t" << yData[n].yz() << "\t"
                << yData[n].zx() << "\t" << yData[n].zy()
                << "\t" << yData[n].zz()
                << endl;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}


void Foam::uspControllerBase::updateTime()
{
    timeData_++;
}


void Foam::uspControllerBase::updateProperties
(
    const dictionary& dict
)
{
    controllerDict_ = dict.subDict("controllerProperties");
    controllerDict_.readIfPresent("controlSwitch", control_);
    controllerDict_.readIfPresent("readStateFromFile", readStateFromFile_);
}


const Foam::word& Foam::uspControllerBase::regionName() const
{
    return regionName_;
}


Foam::scalar Foam::uspControllerBase::density() const
{
    return density_;
}


Foam::scalar& Foam::uspControllerBase::density()
{
    return density_;
}


const Foam::vector& Foam::uspControllerBase::velocity() const
{
    return velocity_;
}


Foam::vector& Foam::uspControllerBase::velocity()
{
    return velocity_;
}


Foam::scalar Foam::uspControllerBase::temperature() const
{
    return temperature_;
}


Foam::scalar& Foam::uspControllerBase::temperature()
{
    return temperature_;
}


Foam::scalar Foam::uspControllerBase::pressure() const
{
    return pressure_;
}


Foam::scalar& Foam::uspControllerBase::pressure()
{
    return pressure_;
}


const Foam::tensor& Foam::uspControllerBase::strainRate() const
{
    return strainRate_;
}


Foam::tensor& Foam::uspControllerBase::strainRate()
{
    return strainRate_;
}


const Foam::vector& Foam::uspControllerBase::tempGradient() const
{
    return tempGradient_;
}


Foam::vector& Foam::uspControllerBase::tempGradient()
{
    return tempGradient_;
}


const Foam::scalarField& Foam::uspControllerBase::densityField() const
{
    return densities_;
}


Foam::scalarField& Foam::uspControllerBase::densityField()
{
    return densities_;
}


const Foam::vectorField& Foam::uspControllerBase::velocityField() const
{
    return velocities_;
}


Foam::vectorField& Foam::uspControllerBase::velocityField()
{
    return velocities_;
}


const Foam::scalarField& Foam::uspControllerBase::temperatureField() const
{
    return temperatures_;
}


Foam::scalarField& Foam::uspControllerBase::temperatureField()
{
    return temperatures_;
}


const Foam::scalarField& Foam::uspControllerBase::pressureField() const
{
    return pressures_;
}


Foam::scalarField& Foam::uspControllerBase::pressureField()
{
    return pressures_;
}


bool Foam::uspControllerBase::singleValueController() const
{
    return singleValueController_;
}


bool& Foam::uspControllerBase::singleValueController()
{
    return singleValueController_;
}


bool Foam::uspControllerBase::fieldController() const
{
    return fieldController_;
}


bool& Foam::uspControllerBase::fieldController()
{
    return fieldController_;
}


bool Foam::uspControllerBase::writeInTimeDir() const
{
    return writeInTimeDir_;
}


bool Foam::uspControllerBase::writeInCase() const
{
    return writeInCase_;
}


Foam::scalar Foam::uspControllerBase::avReqDensity()
{
    return averageProperty(density_, densities_);
}


Foam::vector Foam::uspControllerBase::avReqVelocity()
{
    return averageProperty(velocity_, velocities_);
}


Foam::scalar Foam::uspControllerBase::avReqTemperature()
{
    return averageProperty(temperature_, temperatures_);
}


Foam::scalar Foam::uspControllerBase::avReqPressure()
{
    return averageProperty(pressure_, pressures_);
}


// ************************************************************************* //
