/*---------------------------------------------------------------------------* \
  =========                 |
  \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration     |
    \  /    A nd           | www.openfoam.com
     \/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "uniGasParcel.H"
#include "uniGasCloud.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::uniGasParcel::move
(
    uniGasCloud& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = cloud.pMesh();

    if (newParcel() == 1)
    {
        label& nP = newParcel();
        Random& rndGen = cloud.rndGen();
        stepFraction() = rndGen.sample01<scalar>();
        nP = 0;
    }

    // For reduced - D cases, the velocity used to track needs to be
    // constrained, but the actual U_ of the parcel must not be
    // altered or used, as it is altered by patch interactions an
    // needs to retain its 3D value for collision purposes.
    vector Utracking = U_;

    while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
    {

        Utracking = U_;

        // Apply correction to velocity to constrain tracking for
        // reduced - D cases
        meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);

        // Deviation from the mesh centre for reduced - D cases
        const vector d = deviationFromMeshCentre();

        const scalar f = 1 - stepFraction();
        trackToAndHitFace(f*trackTime*Utracking - d, f, cloud, td);

        if (onFace())
        {
            // apply non-uniform cell weighting
            //if (cloud.cellWeighted())
            //{
            //    cloud.cellWeighting(*this, td.keepParticle);
            //}
            
            //-monitoring flux properties
            cloud.tracker().updateFields(*this);
            cloud.functions().postFace(*this, td.keepParticle);
        }

        if (onBoundaryFace() && td.keepParticle)
        {
            forAll(cloud.boundaries().cyclicBoundaryModels(), c)
            {
                const labelList& faces =
                    cloud.boundaries().cyclicBoundaryModels()[c]->allFaces();

                if (faces.find(this->face()) != -1)
                {
                    cloud.boundaries().cyclicBoundaryModels()[c]->controlMol
                    (
                        *this, td
                    );
                }
            }
        }
    }

    return td.keepParticle;
}


bool Foam::uniGasParcel::hitPatch
(
    uniGasCloud& cloud,
    trackingData& td
)
{
    return false;
}


void Foam::uniGasParcel::hitProcessorPatch
(
    uniGasCloud& cloud,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::uniGasParcel::hitWallPatch
(
    uniGasCloud& cloud,
    trackingData& td
)
{
    const label patchModelI = cloud.boundaries().patchToModelIds()[patch()];

    // Apply a boundary model when a molecule collides with this poly patch
    cloud.boundaries().patchBoundaryModels()[patchModelI]->controlParticle
    (
        *this,
        td
    );
}


void Foam::uniGasParcel::transformProperties
(
    const tensor& T
)
{
   particle::transformProperties(T);
   U_ = transform(T, U_);
}


void Foam::uniGasParcel::transformProperties
(
    const vector& separation
)
{
    particle::transformProperties(separation);
}


// ************************************************************************* //
