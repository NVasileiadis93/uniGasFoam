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

#include "uniGasForceSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasForceSurface, 0);
addToRunTimeSelectionTable(uniGasField, uniGasForceSurface, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasForceSurface::uniGasForceSurface
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
    patchName_(propsDict_.get<word>("patch")),
    averagingAcrossManyRuns_(propsDict_.getOrDefault<bool>("averagingAcrossManyRuns",false)),
    sampleCounter_(0),
    patchId_(-1),
    timeAvCounter_(0.0),
    timeIndex_(0),
    force_(vector::zero),
    forcePatch_(1, vector::zero)
    
{

    // Read stored data from dictionary
    if (averagingAcrossManyRuns_)
    {
        Info << nl << "Averaging across many runs initiated." << nl << endl;
        readIn();
    }

    // Find patchId
    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

    if (patchId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch: " << patchName_ << nl << "in: "
            << mesh_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::uniGasForceSurface::readIn()
{
    localIOdictionary dict
    (
        IOobject
        (
            "forceSurface_" + fieldName_ + "_" + patchName_,
            mesh_.time().timeName(),
            "uniform",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    dict.readIfPresent("force", force_);

}


void Foam::uniGasForceSurface::writeOut()
{
    if (mesh_.time().writeTime())
    {
        localIOdictionary dict
        (
            IOobject
            (
                "forceSurface_" + fieldName_ + "_" + patchName_,
                mesh_.time().timeName(),
                "uniform",
                mesh_.time(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("force", force_);

        dict.regIOobject::writeObject
        (
            IOstreamOption::ASCII,
            dict.time().writeCompression()
        );
    }
}


void Foam::uniGasForceSurface::createField()
{

    Info << "Initialising uniGasForceSurface field" << endl;

}


void Foam::uniGasForceSurface::calculateField()
{

    const scalar& deltaT = mesh_.time().deltaTValue();

    sampleCounter_++;

    if (sampleInterval_ <= sampleCounter_)
    {

        timeAvCounter_ += deltaT;

        // Sample force
        auto& bm = cloud_.boundaryFluxMeasurements();

        forAll(bm.rhoNBF(), i)
        {

            const label iD = typeIds_.find(i);

            if (iD != -1)
            {

                forAll(bm.rhoNBF()[i], j)
                {
                    if( j == patchId_)
                    {
                        forAll(bm.rhoNBF()[i][j], k)
                        {
                            const label& face = mesh_.boundary()[j].start() + k;
                            force_ += bm.fDBF()[i][j][k]*mag(mesh_.faceAreas()[face])*deltaT;
                        }
                    }
                }
            }

        }

        sampleCounter_ = 0;
    }

    // Average force measurement
    if (mesh_.time().writeTime())
    {
      
        if (Pstream::parRun())
        {
            reduce(force_, sumOp<vector>());
        }

        if (patchId_ != -1)
        {
            forcePatch_[timeIndex_] = force_/timeAvCounter_;
        }
        else
        {
            forcePatch_[timeIndex_] = vector::zero;
        }

        if (resetFieldsAtOutput_ && mesh_.time().value() < resetFieldsAtOutputUntilTime_+0.5*deltaT)
        {
            timeAvCounter_ = 0;
            force_ = vector::zero;
        }

        if (averagingAcrossManyRuns_)
        {
            writeOut();
        }

        ++timeIndex_;
    }

}


void Foam::uniGasForceSurface::writeField()
{

    if (mesh_.time().writeTime())
    {
        timeIndex_ = 0;

        if (Pstream::master())
        {
            scalarField timeField(1);

            timeField[0] = mesh_.time().timeOutputValue();

            writeTimeData
            (
                casePath_,
                "force_" + patchName_ + "_"+fieldName_ + ".log",
                timeField,
                forcePatch_,
                true
            );

        }
    }
}


void Foam::uniGasForceSurface::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uniGasField::updateProperties(dict);

    propsDict_.readIfPresent("averagingAcrossManyRuns", averagingAcrossManyRuns_);

}
// ************************************************************************* //
