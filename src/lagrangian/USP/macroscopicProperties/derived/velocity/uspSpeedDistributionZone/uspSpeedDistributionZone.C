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

#include "uspSpeedDistributionZone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspSpeedDistributionZone, 0);

addToRunTimeSelectionTable
(
    uspField,
    uspSpeedDistributionZone,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspSpeedDistributionZone::uspSpeedDistributionZone
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
    distr_(binWidth_)
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

void Foam::uspSpeedDistributionZone::createField()
{
    Info << "Initialising uspSpeedDistributionZone field" << endl;
}


void Foam::uspSpeedDistributionZone::calculateField()
{
    timeVel_++;

    if (timeVel_.samplingTime())
    {
        const auto& cellOccupancy = cloud_.cellOccupancy();

        const labelList& cells = mesh_.cellZones()[regionId_];

        forAll(cells, c)
        {
            const label cellI = cells[c];
            const List<uspParcel*>& parcelsInCell = cellOccupancy[cellI];

            forAll(parcelsInCell, pIC)
            {
                uspParcel* p = parcelsInCell[pIC];

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
            Info << "uspSpeedDistributionZone, averaged velocity "
                << UMean_ << endl;
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

        forAll(cells, c)
        {
            const label cellI = cells[c];

            for (uspParcel* p : cellOccupancy[cellI])
            {
                if (typeIds_.find(p->typeId()) != -1)
                {
                    distr_.add(mag(p->U() - UMean_));
                }
            }
        }
    }
}


void Foam::uspSpeedDistributionZone::writeField()
{
    const Time& runTime = mesh_.time();

    if ((runTime.writeTime()) && (timeVel_.averagingTime()))
    {
        fileName timePath(runTime.path()/runTime.timeName()/"uniform");

        if (!isDir(timePath))
        {
            mkDir(timePath);
        }

        List<Pair<scalar>> normalisedDistriubtion = distr_.normalised();

        label nSize = normalisedDistriubtion.size();

        if (Pstream::parRun())
        {
            reduce(nSize, sumOp<label>());
        }

        scalarField xAxis (nSize, 0.0);
        scalarField yAxis (nSize, 0.0);

        forAll(normalisedDistriubtion, i)
        {
            xAxis[i] = normalisedDistriubtion[i].first();
            yAxis[i] = normalisedDistriubtion[i].second();
        }

        if (Pstream::parRun())
        {
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if (p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::commsTypes::blocking, proc);
                        toNeighbour << xAxis << yAxis;
                    }
                }
            }

            // receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if (p != Pstream::myProcNo())
                {
                    scalarField xAxisProc;
                    scalarField yAxisProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                        fromNeighbour >> xAxisProc >> yAxisProc;
                    }

                    forAll(xAxisProc, i)
                    {
                        xAxis[i] += xAxisProc[i];
                        yAxis[i] += yAxisProc[i];
                    }
                }
            }
        }

        writeTimeData(timePath, "speedDistribution_"+fieldName_+"_"
        +regionName_, xAxis, yAxis);

        if (timeVel_.resetFieldsAtOutput() && (mesh_.time().value() <= timeVel_.resetFieldsAtOutputUntilTime()))
        {
            distr_.clear();
        }
    }
}


void Foam::uspSpeedDistributionZone::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uspField::updateProperties(dict);
}


// ************************************************************************* //
