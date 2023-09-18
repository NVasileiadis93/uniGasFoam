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

#include "uspChapmanEnskogFreeStreamInflowFieldPatch.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspChapmanEnskogFreeStreamInflowFieldPatch, 0);

addToRunTimeSelectionTable
(
    uspGeneralBoundary,
    uspChapmanEnskogFreeStreamInflowFieldPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspChapmanEnskogFreeStreamInflowFieldPatch::uspChapmanEnskogFreeStreamInflowFieldPatch
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    inletTransT_(),
    inletRotT_(),
    inletVibT_(),
    inletElecT_(),
    inletVelocities_(),
    inletHeatFluxes_(),
    inletStresses_(),
    numberDensities_(),
    boundaryTransT_
    (
        IOobject
        (
            "boundaryTransT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    boundaryRotT_
    (
        IOobject
        (
            "boundaryRotT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    boundaryVibT_
    (
        IOobject
        (
            "boundaryVibT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    boundaryElecT_
    (
        IOobject
        (
            "boundaryElecT",
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
    boundaryHeatFlux_
    (
        IOobject
        (
            "boundaryHeatFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    boundaryStress_
    (
        IOobject
        (
            "boundaryStress",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    boundaryNumberDensity_()
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspChapmanEnskogFreeStreamInflowFieldPatch::initialConfiguration()
{}


void Foam::uspChapmanEnskogFreeStreamInflowFieldPatch::calculateProperties()
{}


void Foam::uspChapmanEnskogFreeStreamInflowFieldPatch::controlParcelsBeforeMove()
{
    computeParcelsToInsert
    (
        inletTransT_,
        inletVelocities_,
        inletHeatFluxes_,
        inletStresses_,
        numberDensities_
    );

    insertParcels
    (
        inletTransT_,
        inletRotT_,
        inletVibT_,
        inletElecT_,
        inletVelocities_,
        inletHeatFluxes_,
        inletStresses_,
        numberDensities_
    );
}


void Foam::uspChapmanEnskogFreeStreamInflowFieldPatch::controlParcelsBeforeCollisions()
{}


void Foam::uspChapmanEnskogFreeStreamInflowFieldPatch::controlParcelsAfterCollisions()
{}


void Foam::uspChapmanEnskogFreeStreamInflowFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uspChapmanEnskogFreeStreamInflowFieldPatch::updateProperties
(
    const dictionary& dict
)
{
    // The main properties should be updated first
    uspGeneralBoundary::updateProperties(dict);
}


void Foam::uspChapmanEnskogFreeStreamInflowFieldPatch::setProperties()
{
    // Read in the type ids
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // Set the accumulator

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(), 0.0);
    }

    inletTransT_.setSize(faces_.size());
    inletRotT_.setSize(faces_.size());
    inletVibT_.setSize(faces_.size());
    inletElecT_.setSize(faces_.size());

    forAll(inletTransT_, f)
    {
        inletTransT_[f] = boundaryTransT_.boundaryField()[patchId_][f];
        inletRotT_[f] = boundaryRotT_.boundaryField()[patchId_][f];
        inletVibT_[f] = boundaryVibT_.boundaryField()[patchId_][f];
        inletElecT_[f] = boundaryElecT_.boundaryField()[patchId_][f];
    }

    inletVelocities_.setSize(faces_.size(), Zero);

    forAll(inletVelocities_, f)
    {
        inletVelocities_[f] = boundaryU_.boundaryField()[patchId_][f];
    }

    inletHeatFluxes_.setSize(faces_.size(), Zero);

    forAll(inletHeatFluxes_, f)
    {
        inletHeatFluxes_[f] = boundaryHeatFlux_.boundaryField()[patchId_][f];
    }

    inletStresses_.setSize(faces_.size(), Zero);

    forAll(inletStresses_, f)
    {
        inletStresses_[f] = boundaryStress_.boundaryField()[patchId_][f];
    }

    boundaryNumberDensity_.setSize(typeIds_.size());

    forAll(boundaryNumberDensity_, i)
    {
        const word& moleculeName = cloud_.typeIdList()[typeIds_[i]];

        word nameBoundaryDensity("boundaryNumberDensity_" + moleculeName);

        boundaryNumberDensity_[i].reset
        (
            new volScalarField
            (
                IOobject
                (
                    nameBoundaryDensity,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }

    numberDensities_.setSize(typeIds_.size());

    forAll(numberDensities_, i)
    {
        numberDensities_[i].setSize(faces_.size(), 0.0);

        forAll(numberDensities_[i], f)
        {
            numberDensities_[i][f] =
                boundaryNumberDensity_[i]->boundaryField()[patchId_][f];
        }
    }
}


// ************************************************************************* //
