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

#include "uspMeshFieldFill.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspMeshFieldFill, 0);

addToRunTimeSelectionTable
(
    uspConfiguration,
    uspMeshFieldFill,
    dictionary
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspMeshFieldFill::uspMeshFieldFill(uspCloud& cloud, const dictionary& dict)
:
    uspConfiguration(cloud, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::uspMeshFieldFill::setInitialConfiguration()
{

    Info<< nl << "Initialising particles" << endl;

    //- read initial fields
    volScalarField initialTransT_
    (
        IOobject
        (
            "transT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volScalarField initialRotT_
    (
        IOobject
        (
            "rotT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volScalarField initialVibT_
    (
        IOobject
        (
            "vibT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volScalarField initialElecT_
    (
        IOobject
        (
            "elecT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volVectorField initialU_
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    // Read in the type ids
    List<word> molecules = uspInitialiseDict_.lookup("typeIdList");

    //- List of inlet densities (one entry for each species)
    List<autoPtr<volScalarField>> initialNumberDensityPtr_;

    initialNumberDensityPtr_.setSize(molecules.size());

    forAll(initialNumberDensityPtr_, i)
    {

        const word& moleculeName = molecules[i];

        initialNumberDensityPtr_[i].reset
        (
            new volScalarField
            (
                IOobject
                (
                    "numberDensity_"+moleculeName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

    }

    const auto& meshCC = cloud_.mesh().cellCentres();
    const auto& meshV = cloud_.mesh().V();

    // Initialise subcells in case of adaptive simulation
    if (cloud_.adaptive())
    {
        cloud_.dynamicAdapter().setInitialConfiguration
        (
            initialNumberDensityPtr_,
            initialTransT_
        );
    }

    //Compute cell weights in case of cell weighted simulation    
    if (cloud_.cellWeighted())
    {
        forAll(mesh_.cells(), cell)
        {

            scalar totalNumberDensity = 0.0;
            forAll(molecules, i)
            {
                const volScalarField& initialNumberDensity = initialNumberDensityPtr_[i];
                totalNumberDensity += initialNumberDensity[cell];
            }

            scalar RWF = cloud_.axiRWF(meshCC[cell]);
            const vector& subcellLevels = cloud_.subcellLevels()[cell];
            const scalar nSubcells = subcellLevels.x()*subcellLevels.y()*subcellLevels.z();
            cloud_.cellWeightFactor()[cell] = (totalNumberDensity*meshV[cell])/(cloud_.particlesPerSubcell()*nSubcells*cloud_.nParticle()*RWF);

        }
        cloud_.cellWeightFactor().correctBoundaryConditions();

        scalar maxCellWeightRatio;
        label smoothingPasses = 0;
        do
        {

            smoothingPasses++;

            cloud_.cellWeightFactor() = fvc::average(fvc::interpolate(cloud_.cellWeightFactor()));
            cloud_.cellWeightFactor().correctBoundaryConditions(); 

            forAll(mesh_.cells(), cell)
            {

                scalar totalNumberDensity = 0.0;
                forAll(molecules, i)
                {
                    const volScalarField& initialNumberDensity = initialNumberDensityPtr_[i];
                    totalNumberDensity += initialNumberDensity[cell];
                }

                scalar RWF = cloud_.axiRWF(meshCC[cell]);
                const vector& subcellLevels = cloud_.subcellLevels()[cell];
                const scalar nSubcells = subcellLevels.x()*subcellLevels.y()*subcellLevels.z();
                cloud_.cellWeightFactor()[cell] = 
                    max(min((totalNumberDensity*meshV[cell])/(cloud_.minParticlesPerSubcell()*nSubcells*cloud_.nParticle()*RWF),cloud_.cellWeightFactor()[cell]),SMALL);

            }
            cloud_.cellWeightFactor().correctBoundaryConditions();

            maxCellWeightRatio = VSMALL;
            forAll(mesh_.faces(), face)
            {

                if (mesh_.isInternalFace(face))
                {

                    scalar ownerCWF = cloud_.cellWeightFactor()[mesh_.faceOwner()[face]];
                    scalar neighbourCWF = cloud_.cellWeightFactor()[mesh_.faceNeighbour()[face]];
                    scalar cellWeightRatio = max(ownerCWF/neighbourCWF, neighbourCWF/ownerCWF);
                    if (cellWeightRatio > maxCellWeightRatio)
                    {
                        maxCellWeightRatio = cellWeightRatio;
                    }
                }

            }

            if (Pstream::parRun())
            {
                reduce(maxCellWeightRatio, maxOp<scalar>());
            }

        } while(maxCellWeightRatio > 1.0 + cloud_.maxCellWeightRatio() && smoothingPasses < cloud_.maxSmoothingPasses());
    }

    // Initialise particles
    forAll(mesh_.cells(), cell)
    {
        List<tetIndices> cellTets =
            polyMeshTetDecomposition::cellTetIndices(mesh_, cell);

        for (const tetIndices& cellTetIs : cellTets)
        {
            tetPointRef tet = cellTetIs.tet(mesh_);

            scalar tetVolume = tet.mag();

            forAll(molecules, i)
            {
                const word& moleculeName = molecules[i];

                label typeId(cloud_.typeIdList().find(moleculeName));

                if (typeId == -1)
                {
                    FatalErrorInFunction
                        << "typeId " << moleculeName << "not defined."
                        << abort(FatalError);
                }

                const auto& cP = cloud_.constProps(typeId);

                const volScalarField& initialNumberDensity = initialNumberDensityPtr_[i];

                scalar numberDensity = initialNumberDensity[cell]/cloud_.nParticle();

                scalar translationalTemperature = initialTransT_[cell];

                scalar rotationalTemperature = initialRotT_[cell];

                scalar vibrationalTemperature = initialVibT_[cell];

                scalar electronicTemperature = initialElecT_[cell];

                vector velocity = initialU_[cell];

                // Calculate the number of particles required
                scalar CWF = cloud_.cellWF(cell);
                scalar RWF = cloud_.axiRWF(meshCC[cell]);
                scalar particlesRequired = numberDensity*tetVolume/(CWF*RWF);

                // Only integer numbers of particles can be inserted
                label nParticlesToInsert = label(particlesRequired);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of particlesRequired
                if
                (
                    (particlesRequired - nParticlesToInsert)
                  > rndGen_.sample01<scalar>()
                )
                {
                    ++nParticlesToInsert;
                }

                for (label pI = 0; pI < nParticlesToInsert; ++pI)
                {
                    point p = tet.randomPoint(rndGen_);

                    vector U = vector::zero;

                    U = cloud_.equipartitionLinearVelocity
                    (
                        translationalTemperature,
                        cP.mass()
                    );

                    scalar ERot = cloud_.equipartitionRotationalEnergy
                    (
                        rotationalTemperature,
                        cP.rotationalDoF()
                    );

                    labelList vibLevel =
                        cloud_.equipartitionVibrationalEnergyLevel
                        (
                            vibrationalTemperature,
                            cP.vibrationalDoF(),
                            typeId
                        );

                    label ELevel = cloud_.equipartitionElectronicLevel
                    (
                        electronicTemperature,
                        cP.degeneracyList(),
                        cP.electronicEnergyList(),
                        typeId
                    );

                    U += velocity;

                    label newParcel = 0;

                    scalar CWF = cloud_.cellWF(cell);
                    scalar RWF = cloud_.axiRWF(meshCC[cell]);

                    cloud_.addNewParcel
                    (
                        p,
                        U,
                        CWF,
                        RWF,
                        ERot,
                        ELevel,
                        cell,
                        typeId,
                        newParcel,
                        vibLevel
                    );
                }
            }
        }
    }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the maximum most probable thermal speed (Bird,
    // p222 - 223)

    List<scalar> initialNumberDensityPtr__;
    initialNumberDensityPtr__.setSize(molecules.size());

    forAll(mesh_.cells(), cell)
    {

        forAll(molecules, i)
        {
            const volScalarField& initialNumberDensity = initialNumberDensityPtr_[i];
            initialNumberDensityPtr__[i] = initialNumberDensity[cell];
        }
        
        label mostAbundantType(findMax(initialNumberDensityPtr__));

        const auto& cP = cloud_.constProps(mostAbundantType);

        cloud_.sigmaTcRMax()[cell] =
        cP.sigmaT()
        *cloud_.maxwellianMostProbableSpeed
            (
                initialTransT_[cell],
                cP.mass()
            );

    }

    cloud_.sigmaTcRMax().correctBoundaryConditions();

}


// ************************************************************************* //
