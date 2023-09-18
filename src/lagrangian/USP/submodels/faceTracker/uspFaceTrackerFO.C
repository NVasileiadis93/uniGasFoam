/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "uspFaceTrackerFO.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "uspCloud.H"


namespace Foam
{
makeCloudFunctionObjectTypeNT(uspFaceTrackerFO, uspCloud);
    //defineTypeNameAndDebug(uspFaceTrackerFO, 0);
    //makeCloudFunctionObjectTypeNT(uspFaceTrackerFO, uspCloud);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspFaceTrackerFO::uspFaceTrackerFO
(
    const dictionary& dict,
    uspCloud& owner,
    const word& modelName
)
:
    CloudFunctionObject<uspCloud>(dict, owner, modelName, typeName),
    parcelIdFlux_(owner.typeIdList().size()),
    massIdFlux_(owner.typeIdList().size())
{
    forAll(parcelIdFlux_, i)
    {
        parcelIdFlux_[i].setSize(owner.mesh().nFaces(), 0.0);
        massIdFlux_[i].setSize(owner.mesh().nFaces(), 0.0);
    }
}


Foam::uspFaceTrackerFO::uspFaceTrackerFO(const uspFaceTrackerFO& ft)
:
    CloudFunctionObject<uspCloud>(ft),
    parcelIdFlux_(ft.parcelIdFlux_),
    massIdFlux_(ft.massIdFlux_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspFaceTrackerFO::postEvolve()
{
    forAll(parcelIdFlux_, i)
    {
        parcelIdFlux_[i] = scalar(0);
        massIdFlux_[i] = scalar(0);
    }
}


void Foam::uspFaceTrackerFO::postFace(const uspParcel& p, bool&)
{
    const label crossedFace = p.face();
    const label typeId = p.typeId();
    const auto& constProp = owner().constProps(typeId);
    const scalar mass = constProp.mass();
    const vector& U = p.U();
//     const vector mom = p.U()*mass;

    // Check which patch was hit
    const label patchId = owner().mesh().boundaryMesh().whichPatch(crossedFace);

    // Direction of uspParcel trajectory with respect to the face normal
    scalar sgn = sign(U & owner().mesh().faceAreas()[crossedFace]);

    if (patchId != -1)
    {
        // Boundary face
        const polyPatch& pp = owner().mesh().boundaryMesh()[patchId];

        const label faceIndex = crossedFace - pp.start();

        // Correct cyclic patches
        if (isA<cyclicPolyPatch>(pp))
        {
            label coupledFace =
                refCast<const cyclicPolyPatch>(pp).neighbPatch().start()
              + faceIndex;

            parcelIdFlux_[typeId][coupledFace] += 1.0;
            massIdFlux_[typeId][coupledFace] += mass;
        }

        // Properties are appended to the face of the leaving
        // processor only. Normal vector points out from the domain.
        if (isA<processorPolyPatch>(pp))
        {
            parcelIdFlux_[typeId][crossedFace] += sgn*1.0;
            massIdFlux_[typeId][crossedFace] += sgn*mass;
        }
    }
    else
    {
        // Internal face

        // properties
        parcelIdFlux_[typeId][crossedFace] += sgn*1.0;
        massIdFlux_[typeId][crossedFace] += sgn*mass;
    }
}


// ************************************************************************* //
