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

#include "uspFluxController.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspFluxController, 0);
defineRunTimeSelectionTable(uspFluxController, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspFluxController::uspFluxController
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspControllerBase(mesh, cloud, dict)
{
    const faceZoneMesh& faceZones = mesh_.faceZones();
    regionId_ = faceZones.findZoneID(regionName_);

    if (regionId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find region (faceZone): " << regionName_ << nl << "in: "
            << mesh_.time().system()/"controllersDict"
            << exit(FatalError);
    }

    setFacesInfo();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::uspFluxController> Foam::uspFluxController::New
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("fluxControllerModel"));

    Info<< "Selecting fluxController " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "uspFluxController",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<uspFluxController>(cstrIter()(mesh, cloud, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspFluxController::setFacesInfo()
{
    const labelList& faces = controlZone();

    if (Pstream::parRun())
    {
        DynamicList<label> processorFaces;

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

        processorFaces_.transfer(processorFaces);

        label nInternalFaces = faces.size() - processorFaces.size();
        internalFaces_.setSize(nInternalFaces, -1);

        label counter = 0;

        for (const label faceI : faces)
        {
            if (processorFaces.find(faceI) == -1)
            {
                internalFaces_[counter] = faceI;
                ++counter;
            }
        }

        for (const label faceI : internalFaces_)
        {
            zoneSurfaceArea_ += mag(mesh_.faceAreas()[faceI]);
        }

        // faces on a zone located on a processor cut belong to both
        // processors (hence the 0.5)

        for (const label faceI : processorFaces_)
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
                        OPstream toNeighbour
                        (
                            Pstream::commsTypes::blocking,
                            proc
                        );
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
                        IPstream fromNeighbour
                        (
                            Pstream::commsTypes::blocking,
                            proc
                        );
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
}


const Foam::labelList& Foam::uspFluxController::controlZone() const
{
    return mesh_.faceZones()[regionId_];
}


Foam::label Foam::uspFluxController::isFaceOnControlZone
(
    const label facei
) const
{
    return controlZone().find(facei);
}


// ************************************************************************* //
