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

#include "uniGasMassFluxSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasMassFluxSurface, 0);
addToRunTimeSelectionTable(uniGasField, uniGasMassFluxSurface, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasMassFluxSurface::uniGasMassFluxSurface
(
    const Time& t,
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.get<word>("field")),
    typeIds_(cloud_.getTypeIDs(propsDict_)),
    faceZoneName_(propsDict_.get<word>("faceZone")),        
    fluxDirection_(propsDict_.get<vector>("fluxDirection")),
    averagingAcrossManyRuns_(propsDict_.getOrDefault<bool>("averagingAcrossManyRuns",false)),
    sampleCounter_(0),
    regionId_(-1),
    zoneSurfaceArea_(0.0),
    molsZone_(0.0),
    massZone_(0.0),
    momentumZone_(0.0),
    energyZone_(0.0),
    timeAvCounter_(0.0),
    timeIndex_(0),
    molFluxZone_(1, 0.0),
    massFluxZone_(1, 0.0),
    momentumFluxZone_(1, 0.0),
    energyFluxZone_(1, 0.0)
{

    // Read stored data from dictionary
    if (averagingAcrossManyRuns_)
    {
        Info << nl << "Averaging across many runs initiated." << nl << endl;
        readIn();
    }

    // select face zone
    const faceZoneMesh& faceZones = mesh_.faceZones();

    regionId_ = faceZones.findZoneID(faceZoneName_);

    if (regionId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find region: " << faceZoneName_ << nl << "in: "
            << mesh_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

    fluxDirection_.normalise();

    // find total surface area
    const labelList& faces = faceZones[regionId_];

    if (Pstream::parRun())
    {
        DynamicList<label> processorFaces(0);

        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                for (label p = 0; p < patch.size(); ++p)
                {
                    label patchFaceI = p + patch.start();
                    label faceId = faces.find(patchFaceI);

                    if (faceId != -1)
                    {
                        processorFaces.append(patchFaceI);
                    }
                }
            }
        }

        processorFaces.shrink();

        label nInternalFaces = faces.size() - processorFaces.size();

        List<label> internalFaces(nInternalFaces, 0);

        label counter = 0;

        for (const label faceI : faces)
        {
            if (processorFaces.find(faceI) == -1)
            {
                internalFaces[counter] = faceI;
                ++counter;
            }
        }

        for (const label faceI : internalFaces)
        {
            zoneSurfaceArea_ += mag(mesh_.faceAreas()[faceI]);
        }

        for (const label faceI : processorFaces)
        {
            zoneSurfaceArea_ += 0.5*mag(mesh_.faceAreas()[faceI]);
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
                        toNeighbour << zoneSurfaceArea_;
                    }
                }
            }

            // receiving
            for (int p = 0; p < Pstream::nProcs(); ++p)
            {
                if (p != Pstream::myProcNo())
                {
                    scalar zoneSurfaceAreaProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::commsTypes::blocking, proc);
                        fromNeighbour >> zoneSurfaceAreaProc;
                    }

                    zoneSurfaceArea_ += zoneSurfaceAreaProc;
                }
            }
        }
    }
    else
    {
        for (const label faceI : faces)
        {
            zoneSurfaceArea_ += mag(mesh_.faceAreas()[faceI]);
        }
    }

    Info << "zoneSurfaceArea_ = " << zoneSurfaceArea_ << endl;

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasMassFluxSurface::readIn()
{
    localIOdictionary dict
    (
        IOobject
        (
            "massFluxSurface_" + fieldName_ + "_" + faceZoneName_,
            mesh_.time().timeName(),
            "uniform",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    dict.readIfPresent("timeAvCounter", timeAvCounter_);
    dict.readIfPresent("molsZone", molsZone_);
    dict.readIfPresent("massZone", massZone_);
    dict.readIfPresent("momentumZone", momentumZone_);
    dict.readIfPresent("energyZone", energyZone_);
    
}


void Foam::uniGasMassFluxSurface::writeOut()
{
    if (mesh_.time().writeTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "massFluxSurface_"+fieldName_+"_"+faceZoneName_,
                mesh_.time().timeName(),
                "uniform",
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("timeAvCounter", timeAvCounter_);
        dict.add("molsZone", molsZone_);
        dict.add("massZone", massZone_);
        dict.add("momentumZone", momentumZone_);
        dict.add("energyZone", energyZone_);

        dict.regIOobject::writeObject
        (
            IOstreamOption(mesh_.time().writeFormat()),
            true
        );
    }
}


void Foam::uniGasMassFluxSurface::createField()
{

    Info << "Initialising uniGasMassFluxSurface field" << endl;

}


void Foam::uniGasMassFluxSurface::calculateField()
{

const scalar& deltaT = mesh_.time().deltaTValue();

    sampleCounter_++;

    if (sampleInterval_ <= sampleCounter_)
    {

        timeAvCounter_ += deltaT;

        const List<scalarField>& molIdFlux = cloud_.tracker().parcelIdFlux();
        const List<scalarField>& massIdFlux = cloud_.tracker().massIdFlux();
        const List<vectorField>& momentumIdFlux = cloud_.tracker().momentumIdFlux();
        const List<scalarField>& energyIdFlux = cloud_.tracker().energyIdFlux();

        scalar molFlux = 0.0;
        scalar massFlux = 0.0;
        scalar momentumFlux = 0.0;
        scalar energyFlux = 0.0;

        const faceZoneMesh& faceZones = mesh_.faceZones();
        const labelList& faces = faceZones[regionId_];

        for (const label faceI : faces)
        {

            vector nF = mesh_.faceAreas()[faceI]/mag(mesh_.faceAreas()[faceI]);

            forAll(molIdFlux, i)
            {

                const label iD = typeIds_.find(i);

                if (iD != -1)
                {

                    const scalar& nParticle = cloud_.nParticle();

                    molFlux += (molIdFlux[i][faceI]*nParticle*nF) & fluxDirection_;
                    massFlux += (massIdFlux[i][faceI]*nParticle*nF) & fluxDirection_;
                    momentumFlux += (momentumIdFlux[i][faceI]*nParticle) & fluxDirection_;
                    energyFlux += (energyIdFlux[i][faceI]*nParticle*nF) & fluxDirection_;
                }
            }

        }

        molsZone_ += molFlux;
        massZone_ += massFlux;
        momentumZone_ += momentumFlux;
        energyZone_ += energyFlux;

        sampleCounter_ = 0;
    }

    const Time& runTime = mesh_.time();

    // Average measurement and calculate properties
    if (runTime.writeTime())
    {
        
        scalar molsZone = molsZone_;
        scalar massZone = massZone_;
        scalar momentumZone = momentumZone_;
        scalar energyZone = energyZone_;

        if (Pstream::parRun())
        {
            reduce(molsZone, sumOp<scalar>());
            reduce(massZone, sumOp<scalar>());
            reduce(momentumZone, sumOp<scalar>());
            reduce(energyZone, sumOp<scalar>());
        }

        if (zoneSurfaceArea_ > 0.0)
        {
            molFluxZone_[timeIndex_] = molsZone/(timeAvCounter_*zoneSurfaceArea_);
            massFluxZone_[timeIndex_] = massZone/(timeAvCounter_*zoneSurfaceArea_);
            momentumFluxZone_[timeIndex_] = momentumZone/(timeAvCounter_*zoneSurfaceArea_);
            energyFluxZone_[timeIndex_] = energyZone/(timeAvCounter_*zoneSurfaceArea_);
        }
        else
        {
            molFluxZone_[timeIndex_] = 0.0;
            massFluxZone_[timeIndex_] = 0.0;
            momentumFluxZone_[timeIndex_] = 0.0;
            energyFluxZone_[timeIndex_] = 0.0;
        }

        if (resetFieldsAtOutput_ && mesh_.time().value() < resetFieldsAtOutputUntilTime_+0.5*deltaT)
        {
            timeAvCounter_ = 0;
            molsZone_ = 0.0;
            massZone_ = 0.0;
            momentumZone_ = 0.0;
            energyZone_ = 0.0;
        }

        if (averagingAcrossManyRuns_)
        {
            writeOut();
        }

        ++timeIndex_;
    }

}


void Foam::uniGasMassFluxSurface::writeField()
{
    const Time& runTime = mesh_.time();

    if (runTime.writeTime())
    {
        timeIndex_ = 0;

        if (Pstream::master())
        {
            scalarField timeField(1);

            timeField[0] = mesh_.time().timeOutputValue();

            writeTimeData
            (
                casePath_,
                "particleFlux_" + faceZoneName_ + "_"+fieldName_ + ".log",
                timeField,
                molFluxZone_,
                true
            );

            writeTimeData
            (
                casePath_,
                "massFlux_" + faceZoneName_ + "_" + fieldName_ + ".log",
                timeField,
                massFluxZone_,
                true
            );

            writeTimeData
            (
                casePath_,
                "momentumFlux_" + faceZoneName_ + "_"+fieldName_ + ".log",
                timeField,
                momentumFluxZone_,
                true
            );

            writeTimeData
            (
                casePath_,
                "energyFlux_" + faceZoneName_ + "_"+fieldName_ + ".log",
                timeField,
                energyFluxZone_,
                true
            );

        }
    }
}


void Foam::uniGasMassFluxSurface::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uniGasField::updateProperties(dict);

    propsDict_.readIfPresent("averagingAcrossManyRuns", averagingAcrossManyRuns_);

    if (averagingAcrossManyRuns_)
    {
        Info << "averagingAcrossManyRuns initiated." << endl;
    }

}
// ************************************************************************* //
