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

#include "uspSpecularWallPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspSpecularWallPatch, 0);

addToRunTimeSelectionTable
(
    uspPatchBoundary,
    uspSpecularWallPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspSpecularWallPatch::uspSpecularWallPatch
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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspSpecularWallPatch::initialConfiguration()
{}


void Foam::uspSpecularWallPatch::calculateProperties()
{}


void Foam::uspSpecularWallPatch::controlParticle
(
    uspParcel& p,
    uspParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    specularReflection(p);

    measurePropertiesAfterControl(p, 0.0);
}


void Foam::uspSpecularWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uspSpecularWallPatch::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uspPatchBoundary::updateProperties(dict);
}


// ************************************************************************* //
