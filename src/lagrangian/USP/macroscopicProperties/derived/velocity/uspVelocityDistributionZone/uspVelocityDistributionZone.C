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

#include "uspVelocityDistributionZone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspVelocityDistributionZone, 0);

addToRunTimeSelectionTable
(
    uspField,
    uspVelocityDistributionZone,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspVelocityDistributionZone::uspVelocityDistributionZone
(
    const Time& t,
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.get<word>("field")),
    typeIds_(),
    timeDictVel_(dict.subDict("timePropertiesForVelocity")),
    timeVel_(t, timeDictVel_),
    regionName_(propsDict_.get<word>("zone")),
    regionId_(-1),
    UMean_(Zero),
    Ucollected_(Zero),
    nParcels_(0),
    binWidth_(propsDict_.get<scalar>("binWidth")),
    distrX_(binWidth_),
    distrY_(binWidth_),
    distrZ_(binWidth_)
{
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if (regionId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find region: " << regionName_ << nl << "in: "
            << mesh_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspVelocityDistributionZone::createField()
{
    Info<< "Initialising uspVelocityDistributionZone field" << endl;
}


void Foam::uspVelocityDistributionZone::calculateField()
{
    timeVel_++;

    if (timeVel_.samplingTime())
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();

        const labelList& cells = mesh_.cellZones()[regionId_];

        forAll(cells, c)
        {
            const label cellI = cells[c];

            for (uspParcel* p : cellOccupancy[cellI])
            {
                if (typeIds_.find(p->typeId()) != -1)
                {
                    ++nParcels_;
                    Ucollected_ += p->U();
                }
            }
        }
    }

    if (timeVel_.averagingTime())
    {
        vector Ucollected = Ucollected_;
        label nParcels = nParcels_;

        if (Pstream::parRun())
        {
            reduce(Ucollected, sumOp<vector>());
            reduce(nParcels, sumOp<label>());
        }

        if (nParcels_ > 0)
        {
            UMean_ = Ucollected/nParcels;
            Info<< "uspVelocityDistributionZone, averaged velocity " << UMean_ << endl;
        }

        if (timeVel_.resetFieldsAtOutput() && (mesh_.time().value() <= timeVel_.resetFieldsAtOutputUntilTime()))
        {
            Ucollected_ = vector::zero;
            nParcels_ = 0;
        }
    }

    if (timeVel_.samplingTime())
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();

        const labelList& cells = mesh_.cellZones()[regionId_];

        for (const label cellI : cells)
        {
            for (uspParcel* p : cellOccupancy[cellI])
            {
                if (typeIds_.find(p->typeId()) != -1)
                {
                    distrX_.add((p->U().x() - UMean_.x()));
                    distrY_.add((p->U().y() - UMean_.y()));
                    distrZ_.add((p->U().z() - UMean_.z()));
                }
            }
        }
    }
}


void Foam::uspVelocityDistributionZone::writeField()
{
    const Time& runTime = mesh_.time();

    if ((runTime.writeTime()) && (timeVel_.averagingTime()))
    {
        fileName timePath(runTime.path()/runTime.timeName()/"uniform");

        if (!isDir(timePath))
        {
            mkDir(timePath);
        }

        List<Pair<scalar>> rawDistriubtionX(distrX_.raw());
        List<Pair<scalar>> rawDistriubtionY(distrY_.raw());
        List<Pair<scalar>> rawDistriubtionZ(distrZ_.raw());

        label nSizeX = rawDistriubtionX.size();
        label nSizeY = rawDistriubtionY.size();
        label nSizeZ = rawDistriubtionZ.size();

        if (Pstream::parRun())
        {
            reduce(nSizeX, sumOp<label>());
            reduce(nSizeY, sumOp<label>());
            reduce(nSizeZ, sumOp<label>());
        }

        scalarField xAxisX(nSizeX, Zero);
        scalarField yAxisX(nSizeX, Zero);
        scalarField xAxisY(nSizeY, Zero);
        scalarField yAxisY(nSizeY, Zero);
        scalarField xAxisZ(nSizeZ, Zero);
        scalarField yAxisZ(nSizeZ, Zero);


        forAll(rawDistriubtionX, i)
        {
            xAxisX[i] = rawDistriubtionX[i].first();
            yAxisX[i] = rawDistriubtionX[i].second();
        }

        forAll(rawDistriubtionY, i)
        {
            xAxisY[i] = rawDistriubtionY[i].first();
            yAxisY[i] = rawDistriubtionY[i].second();
        }

        forAll(rawDistriubtionZ, i)
        {
            xAxisZ[i] = rawDistriubtionZ[i].first();
            yAxisZ[i] = rawDistriubtionZ[i].second();
        }

        if (Pstream::parRun())
        {
            for (int p = 0; p < Pstream::nProcs(); ++p)
            {
                if (p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                        toNeighbour
                            << xAxisX << yAxisX << xAxisY
                            << yAxisY << xAxisZ << yAxisZ;
                    }
                }
            }

            // receiving
            for (int p = 0; p < Pstream::nProcs(); ++p)
            {
                if (p != Pstream::myProcNo())
                {
                    scalarField xAxisXProc;
                    scalarField yAxisXProc;
                    scalarField xAxisYProc;
                    scalarField yAxisYProc;
                    scalarField xAxisZProc;
                    scalarField yAxisZProc;

                    {
                        IPstream fromNeighbour(Pstream::commsTypes::blocking, p);
                        fromNeighbour
                            >> xAxisXProc >> yAxisXProc >> xAxisYProc
                            >> yAxisYProc >> xAxisZProc >> yAxisZProc;
                    }

                    forAll(xAxisXProc, i)
                    {
                        xAxisX[i] += xAxisXProc[i];
                        yAxisX[i] += yAxisXProc[i];
                    }

                    forAll(xAxisYProc, i)
                    {
                        xAxisY[i] += xAxisYProc[i];
                        yAxisY[i] += yAxisYProc[i];
                    }

                    forAll(xAxisZProc, i)
                    {
                        xAxisZ[i] += xAxisZProc[i];
                        yAxisZ[i] += yAxisZProc[i];
                    }
                }
            }
        }

        const word suffix(fieldName_ + "_" + regionName_);
        writeTimeData
        (
            timePath,
            "velocityDistributionX_" + suffix,
            xAxisX,
            yAxisX
        );
        writeTimeData
        (
            timePath,
            "velocityDistributionY_" + suffix,
            xAxisY,
            yAxisY
        );
        writeTimeData
        (
            timePath,
            "velocityDistributionZ_" + suffix,
            xAxisZ,
            yAxisZ
        );

        if (timeVel_.resetFieldsAtOutput() && (mesh_.time().value() <= timeVel_.resetFieldsAtOutputUntilTime()))
        {
            distrX_.clear();
            distrY_.clear();
            distrZ_.clear();
        }
    }
}


void Foam::uspVelocityDistributionZone::updateProperties
(
    const dictionary& dict
)
{
    // The main properties should be updated first
    uspField::updateProperties(dict);
}


// ************************************************************************* //
