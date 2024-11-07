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

#include "uniGasBoundaryBase.H"
#include "uniGasCloud.H"
#include "uniGasBoundaries.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasBoundaryBase, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasBoundaryBase::uniGasBoundaryBase
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict,
    const word& patchType
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    time_(mesh_.time()),
    patchType_(patchType),
    boundaryDict_(dict.subDict(patchType & "BoundaryProperties")),

    patchName_(boundaryDict_.get<word>("patch")),
    patchId_(mesh_.boundaryMesh().findPatchID(patchName_)),

    numberDensities_(),
    densities_(),
    velocities_(),
    temperatures_(),
    writeInTimeDir_(true),
    writeInCase_(true)
{
    // Confirm that the patch exists on the mesh
    if (patchId_ == -1)
    {
        FatalIOErrorInFunction(boundaryDict_)
            << "Cannot find patch: " << patchName_ << nl << "in: "
            << mesh_.time().system()/uniGasBoundaries::dictName
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasBoundaryBase::updateProperties(const dictionary& dict)
{
    boundaryDict_ = dict.subDict(patchType_ & "BoundaryProperties");
}


void Foam::uniGasBoundaryBase::writeTimeData
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


void Foam::uniGasBoundaryBase::writeTimeData
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


const Foam::word& Foam::uniGasBoundaryBase::patchName() const
{
    return patchName_;
}


Foam::label Foam::uniGasBoundaryBase::patchId() const
{
    return patchId_;
}


const Foam::scalarField& Foam::uniGasBoundaryBase::numberDensities() const
{
    return numberDensities_;
}


Foam::scalarField& Foam::uniGasBoundaryBase::numberDensities()
{
    return numberDensities_;
}


const Foam::scalarField& Foam::uniGasBoundaryBase::densityField() const
{
    return densities_;
}


Foam::scalarField& Foam::uniGasBoundaryBase::densityField()
{
    return densities_;
}


const Foam::vectorField& Foam::uniGasBoundaryBase::velocityField() const
{
    return velocities_;
}


Foam::vectorField& Foam::uniGasBoundaryBase::velocityField()
{
    return velocities_;
}


const Foam::scalarField& Foam::uniGasBoundaryBase::temperatureField() const
{
    return temperatures_;
}


Foam::scalarField& Foam::uniGasBoundaryBase::temperatureField()
{
    return temperatures_;
}


const Foam::tensor& Foam::uniGasBoundaryBase::strainRate() const
{
    return strainRate_;
}


Foam::tensor& Foam::uniGasBoundaryBase::strainRate()
{
    return strainRate_;
}


bool Foam::uniGasBoundaryBase::writeInTimeDir() const
{
    return writeInTimeDir_;
}


bool Foam::uniGasBoundaryBase::writeInCase() const
{
    return writeInCase_;
}


// ************************************************************************* //
