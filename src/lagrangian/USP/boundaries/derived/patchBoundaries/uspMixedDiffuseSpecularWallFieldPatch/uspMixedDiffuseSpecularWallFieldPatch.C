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

#include "uspMixedDiffuseSpecularWallFieldPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspMixedDiffuseSpecularWallFieldPatch, 0);

addToRunTimeSelectionTable
(
    uspPatchBoundary,
    uspMixedDiffuseSpecularWallFieldPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspMixedDiffuseSpecularWallFieldPatch::
uspMixedDiffuseSpecularWallFieldPatch
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    boundaryT_
    (
        IOobject
        (
            "boundaryT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    boundaryU_
    (
        IOobject
        (
            "boundaryU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    diffuseFraction_(propsDict_.get<scalar>("diffuseFraction"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspMixedDiffuseSpecularWallFieldPatch::initialConfiguration()
{}


void Foam::uspMixedDiffuseSpecularWallFieldPatch::calculateProperties()
{}


void Foam::uspMixedDiffuseSpecularWallFieldPatch::controlParticle
(
    uspParcel& p,
    uspParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    label wppIndex = p.patch();
    const polyPatch& patch = mesh_.boundaryMesh()[wppIndex];
    label wppLocalFace = patch.whichFace(p.face());

    if (diffuseFraction_ > cloud_.rndGen().sample01<scalar>())
    {
        // Diffuse reflection
        scalar T = boundaryT_.boundaryField()[wppIndex][wppLocalFace];
        vector U = boundaryU_.boundaryField()[wppIndex][wppLocalFace];

        diffuseReflection(p, T, U);
    }
    else
    {
        // Specular reflection
        specularReflection(p);
    }

    measurePropertiesAfterControl(p, 0.0);
}


void Foam::uspMixedDiffuseSpecularWallFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uspMixedDiffuseSpecularWallFieldPatch::updateProperties
(
    const dictionary& dict
)
{
    // The main properties should be updated first
    uspPatchBoundary::updateProperties(dict);
}


// ************************************************************************* //
