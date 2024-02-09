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

#include "uspDeletionPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspDeletionPatch, 0);
addToRunTimeSelectionTable(uspPatchBoundary, uspDeletionPatch, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspDeletionPatch::uspDeletionPatch
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspPatchBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    allSpecies_(propsDict_.getOrDefault<bool>("allSpecies",false)),
    typeIds_()
{
    measurePropertiesAtWall_ = false;
    writeInTimeDir_ = false;
    writeInCase_ = true;

    if (!allSpecies_)
    {
        typeIds_ = cloud_.getTypeIDs(propsDict_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspDeletionPatch::initialConfiguration()
{}


void Foam::uspDeletionPatch::calculateProperties()
{}


void Foam::uspDeletionPatch::controlParticle
(
    uspParcel& p,
    uspParcel::trackingData& td
)
{
    td.keepParticle = false;
}


void Foam::uspDeletionPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uspDeletionPatch::updateProperties(const dictionary& dict)
{
    // the main properties should be updated first
    uspPatchBoundary::updateProperties(dict);

    propsDict_ = dict.subDict(typeName + "Properties");

    propsDict_.readIfPresent("allSpecies", allSpecies_);

    if (!allSpecies_)
    {
        typeIds_ = cloud_.getTypeIDs(propsDict_);
    }

}


// ************************************************************************* //
