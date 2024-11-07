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

#include "uniGasMeshFill.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasMeshFill, 0);

addToRunTimeSelectionTable
(
    uniGasConfiguration,
    uniGasMeshFill,
    dictionary
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasMeshFill::uniGasMeshFill(uniGasCloud& cloud, const dictionary& dict)
:
    uniGasConfiguration(cloud, dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasMeshFill::setInitialConfiguration()
{
    Info<< nl << "Initialising particles" << endl;

    const scalar translationalTemperature
    (
        uniGasInitialisationDict_.get<scalar>("translationalTemperature")
    );

    const scalar rotationalTemperature
    (
        uniGasInitialisationDict_.get<scalar>("rotationalTemperature")
    );

    const scalar vibrationalTemperature
    (
        uniGasInitialisationDict_.get<scalar>("vibrationalTemperature")
    );

    const scalar electronicTemperature
    (
        uniGasInitialisationDict_.get<scalar>("electronicTemperature")
    );

    const vector velocity(uniGasInitialisationDict_.get<vector>("velocity"));

    const dictionary& numberDensitiesDict = uniGasInitialisationDict_.subDict("numberDensities");

    List<word> molecules(numberDensitiesDict.toc());

    Field<scalar> numberDensities(molecules.size());

    scalar totalNumberDensity  = 0.0;
    forAll(molecules, i)
    {
        numberDensities[i] = numberDensitiesDict.get<scalar>(molecules[i]);
        totalNumberDensity += numberDensities[i];
    }

    const auto& meshCC = cloud_.mesh().cellCentres();
    const auto& meshV = cloud_.mesh().V();

    // Initialise subCells and time-step in case of adaptive simulation
    if (cloud_.adaptive())
    {
        cloud_.dynamicAdapter().setInitialConfiguration
        (
            numberDensities,
            translationalTemperature,
            velocity
        );
    }

    //Compute cell weights in case of cell weighted simulation    
    if (cloud_.cellWeighted())
    {
        forAll(mesh_.cells(), cell)
        {

            scalar RWF = cloud_.axiRWF(meshCC[cell]);
            const vector& subCellLevels = cloud_.subCellLevels()[cell];
            const scalar nSubCells = subCellLevels.x()*subCellLevels.y()*subCellLevels.z();
            cloud_.cellWeightFactor()[cell] = (totalNumberDensity*meshV[cell])/(cloud_.particlesPerSubCell()*nSubCells*cloud_.nParticle()*RWF);

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

                scalar RWF = cloud_.axiRWF(meshCC[cell]);
                const vector& subCellLevels = cloud_.subCellLevels()[cell];
                const scalar nSubCells = subCellLevels.x()*subCellLevels.y()*subCellLevels.z();
                cloud_.cellWeightFactor()[cell] = 
                    max(min((totalNumberDensity*meshV[cell])/(cloud_.minParticlesPerSubCell()*nSubCells*cloud_.nParticle()*RWF),cloud_.cellWeightFactor()[cell]),SMALL);

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

                scalar numberDensity = numberDensities[i]/cloud_.nParticle();

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
    // the most abundant species and the most probable thermal speed (Bird,
    // p222 - 223)

    label mostAbundantType(findMax(numberDensities));

    const auto& cP = cloud_.constProps(mostAbundantType);

    cloud_.sigmaTcRMax().primitiveFieldRef() =
        cP.sigmaT()
       *cloud_.maxwellianMostProbableSpeed
        (
            translationalTemperature,
            cP.mass()
        );

    cloud_.sigmaTcRMax().correctBoundaryConditions();
}


// ************************************************************************* //
