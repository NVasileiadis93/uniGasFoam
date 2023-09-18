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

#include "uspMixedDiffuseSpecularWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspMixedDiffuseSpecularWallPatch, 0);

addToRunTimeSelectionTable
(
    uspPatchBoundary,
    uspMixedDiffuseSpecularWallPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspMixedDiffuseSpecularWallPatch::uspMixedDiffuseSpecularWallPatch
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    diffuseFraction_(propsDict_.get<scalar>("diffuseFraction"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspMixedDiffuseSpecularWallPatch::initialConfiguration()
{}


void Foam::uspMixedDiffuseSpecularWallPatch::calculateProperties()
{}


void Foam::uspMixedDiffuseSpecularWallPatch::controlParticle
(
    uspParcel& p,
    uspParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    if (diffuseFraction_ > cloud_.rndGen().sample01<scalar>())
    {
        // Diffuse reflection
        diffuseReflection(p, temperature_, velocity_);
    }
    else
    {
        // Specular reflection
        specularReflection(p);
    }

    measurePropertiesAfterControl(p, 0.0);
}


void Foam::uspMixedDiffuseSpecularWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uspMixedDiffuseSpecularWallPatch::updateProperties
(
    const dictionary& dict
)
{
    // The main properties should be updated first
    uspPatchBoundary::updateProperties(dict);

    setProperties();
}


void Foam::uspMixedDiffuseSpecularWallPatch::setProperties()
{
    velocity_ = propsDict_.get<vector>("velocity");
    temperature_ = propsDict_.get<scalar>("temperature");
}


// ************************************************************************* //
