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

#include "uspDiffuseWallFieldPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspDiffuseWallFieldPatch, 0);

addToRunTimeSelectionTable
(
    uspPatchBoundary,
    uspDiffuseWallFieldPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspDiffuseWallFieldPatch::uspDiffuseWallFieldPatch
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
    )
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspDiffuseWallFieldPatch::initialConfiguration()
{}


void Foam::uspDiffuseWallFieldPatch::calculateProperties()
{}


void Foam::uspDiffuseWallFieldPatch::controlParticle
(
    uspParcel& p,
    uspParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    label wppIndex = p.patch();
    const polyPatch& patch = mesh_.boundaryMesh()[wppIndex];
    label wppLocalFace = patch.whichFace(p.face());

    scalar T = boundaryT_.boundaryField()[wppIndex][wppLocalFace];
    vector U = boundaryU_.boundaryField()[wppIndex][wppLocalFace];

    diffuseReflection(p, T, U);

    measurePropertiesAfterControl(p, 0.0);
}


void Foam::uspDiffuseWallFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uspDiffuseWallFieldPatch::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uspPatchBoundary::updateProperties(dict);
}


// ************************************************************************* //
