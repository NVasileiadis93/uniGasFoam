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

#include "uspDiffuseWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspDiffuseWallPatch, 0);

addToRunTimeSelectionTable
(
    uspPatchBoundary,
    uspDiffuseWallPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspDiffuseWallPatch::uspDiffuseWallPatch
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspDiffuseWallPatch::initialConfiguration()
{}


void Foam::uspDiffuseWallPatch::calculateProperties()
{}


void Foam::uspDiffuseWallPatch::controlParticle
(
    uspParcel& p,
    uspParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);
    diffuseReflection(p, temperature_, velocity_);
    measurePropertiesAfterControl(p, 0.0);
}


void Foam::uspDiffuseWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uspDiffuseWallPatch::updateProperties(const dictionary& dict)
{
    // the main properties should be updated first
    uspPatchBoundary::updateProperties(dict);

    propsDict_ = dict.subDict(typeName + "Properties");

    setProperties();
}


void Foam::uspDiffuseWallPatch::setProperties()
{
    velocity_ = propsDict_.get<vector>("velocity");
    temperature_ = propsDict_.get<scalar>("temperature");
}


// ************************************************************************* //
