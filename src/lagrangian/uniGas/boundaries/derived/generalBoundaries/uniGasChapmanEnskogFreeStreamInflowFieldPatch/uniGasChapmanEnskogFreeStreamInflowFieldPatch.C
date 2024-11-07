/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "uniGasChapmanEnskogFreeStreamInflowFieldPatch.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasChapmanEnskogFreeStreamInflowFieldPatch, 0);

addToRunTimeSelectionTable
(
    uniGasGeneralBoundary,
    uniGasChapmanEnskogFreeStreamInflowFieldPatch,
    dictionary
);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasChapmanEnskogFreeStreamInflowFieldPatch::uniGasChapmanEnskogFreeStreamInflowFieldPatch
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasGeneralBoundary(mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    numberDensities_(),
    inletTransT_(),
    inletRotT_(),
    inletVibT_(),
    inletElecT_(),
    inletVelocities_(),
    inletHeatFluxes_(),
    inletStresses_(),
    boundaryNumberDensity_(),
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
    )
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    // Get type IDs
    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // Set the accumulator
    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(faces_.size(), 0.0);
    }

    // Set macro properties
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

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasChapmanEnskogFreeStreamInflowFieldPatch::initialConfiguration()
{}


void Foam::uniGasChapmanEnskogFreeStreamInflowFieldPatch::calculateProperties()
{}


void Foam::uniGasChapmanEnskogFreeStreamInflowFieldPatch::controlParcelsBeforeMove()
{
    computeParcelsToInsert
    (
        numberDensities_,
        inletTransT_,
        inletVelocities_,
        inletHeatFluxes_,
        inletStresses_
    );

    insertParcels
    (
        numberDensities_,
        inletTransT_,
        inletRotT_,
        inletVibT_,
        inletElecT_,
        inletVelocities_,
        inletHeatFluxes_,
        inletStresses_
    );
}


void Foam::uniGasChapmanEnskogFreeStreamInflowFieldPatch::controlParcelsBeforeCollisions()
{}


void Foam::uniGasChapmanEnskogFreeStreamInflowFieldPatch::controlParcelsAfterCollisions()
{}


void Foam::uniGasChapmanEnskogFreeStreamInflowFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void Foam::uniGasChapmanEnskogFreeStreamInflowFieldPatch::updateProperties
(
    const dictionary& dict
)
{
    // The main properties should be updated first
    uniGasGeneralBoundary::updateProperties(dict);
}

// ************************************************************************* //
