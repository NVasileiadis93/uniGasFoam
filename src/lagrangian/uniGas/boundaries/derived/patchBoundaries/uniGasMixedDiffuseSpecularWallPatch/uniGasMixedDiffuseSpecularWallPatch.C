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

#include "uniGasMixedDiffuseSpecularWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasMixedDiffuseSpecularWallPatch, 0);

addToRunTimeSelectionTable
(
    uniGasPatchBoundary,
    uniGasMixedDiffuseSpecularWallPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasMixedDiffuseSpecularWallPatch::uniGasMixedDiffuseSpecularWallPatch
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    diffuseFraction_(propsDict_.get<scalar>("diffuseFraction")),
    temperature_(propsDict_.get<scalar>("temperature")),
    velocity_(propsDict_.get<vector>("velocity"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasMixedDiffuseSpecularWallPatch::initialConfiguration()
{}


void Foam::uniGasMixedDiffuseSpecularWallPatch::calculateProperties()
{}


void Foam::uniGasMixedDiffuseSpecularWallPatch::controlParticle
(
    uniGasParcel& p,
    uniGasParcel::trackingData& td
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


void Foam::uniGasMixedDiffuseSpecularWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uniGasMixedDiffuseSpecularWallPatch::updateProperties
(
    const dictionary& dict
)
{
    // The main properties should be updated first
    uniGasPatchBoundary::updateProperties(dict);

    propsDict_ = dict.subDict(typeName + "Properties");

    temperature_ = propsDict_.get<scalar>("temperature");

    velocity_ = propsDict_.get<vector>("velocity");

}

// ************************************************************************* //
