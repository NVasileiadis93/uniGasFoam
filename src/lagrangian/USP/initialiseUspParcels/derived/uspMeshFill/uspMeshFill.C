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

#include "uspMeshFill.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspMeshFill, 0);

addToRunTimeSelectionTable
(
    uspConfiguration,
    uspMeshFill,
    dictionary
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspMeshFill::uspMeshFill(uspCloud& cloud, const dictionary& dict)
:
    uspConfiguration(cloud, dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspMeshFill::setInitialConfiguration()
{
    Info<< nl << "Initialising particles" << endl;

    const scalar translationalTemperature
    (
        uspInitialiseDict_.get<scalar>("translationalTemperature")
    );

    const scalar rotationalTemperature
    (
        uspInitialiseDict_.get<scalar>("rotationalTemperature")
    );

    const scalar vibrationalTemperature
    (
        uspInitialiseDict_.get<scalar>("vibrationalTemperature")
    );

    const scalar electronicTemperature
    (
        uspInitialiseDict_.get<scalar>("electronicTemperature")
    );

    const vector velocity(uspInitialiseDict_.get<vector>("velocity"));

    const dictionary& numberDensitiesDict = uspInitialiseDict_.subDict("numberDensities");

    List<word> molecules(numberDensitiesDict.toc());

    Field<scalar>numberDensities(molecules.size());

    scalar totalNumberDensity  = 0.0;
    forAll(molecules, i)
    {
        numberDensities[i] = numberDensitiesDict.get<scalar>(molecules[i]);
        totalNumberDensity += numberDensities[i];
    }

    const auto& meshCC = cloud_.mesh().cellCentres();
    const auto& meshV = cloud_.mesh().V();

    // Compute subcell levels
    if (cloud_.dynamicAdapter().subcellAdaptation())
    {
        const boolVector& solutionDimensions = cloud_.solutionDimensions(); 
        scalarList speciesMFP(molecules.size(), 0.0);
        scalarList speciesMCT(molecules.size(), 0.0);       

        forAll(mesh_.cells(), cell)
        {

            speciesMFP = 0.0;
            speciesMCT = 0.0; 

            forAll(molecules, i)
            {
                label qspec = 0;

                for (qspec=0; qspec<molecules.size(); ++qspec)
                {
                    scalar dPQ = 0.5*(cloud_.constProps(i).d() + cloud_.constProps(qspec).d());

                    scalar omegaPQ = 0.5*(cloud_.constProps(i).omega() + cloud_.constProps(qspec).omega());

                    scalar massRatio = cloud_.constProps(i).mass()/cloud_.constProps(qspec).mass();

                    if (numberDensities[qspec] > VSMALL)
                    {

                        scalar nDensQ = numberDensities[qspec]/meshV[cell];

                        scalar reducedMass =
                            cloud_.constProps(i).mass()*cloud_.constProps(qspec).mass()
                            /(cloud_.constProps(i).mass() + cloud_.constProps(qspec).mass());

                        //Bird, eq (4.76)
                        speciesMFP[i] += pi*dPQ*dPQ*nDensQ*pow(cloud_.collTref()/translationalTemperature, omegaPQ-0.5)*sqrt(1.0 + massRatio); 

                        // //Bird, eq (4.74)
                        speciesMCT[i] += 
                            2.0*sqrt(pi)*dPQ*dPQ*nDensQ*pow(translationalTemperature/cloud_.collTref(),1.0-omegaPQ)*sqrt(2.0*physicoChemical::k.value()*cloud_.collTref()/reducedMass); 
                    }
                }

            }

            scalar MFP = 0.0;
            scalar MCT = 0.0;
            forAll(molecules, i)
            {
                if (numberDensities[i] > VSMALL)
                {
                    scalar nDensP = numberDensities[i]/meshV[cell];

                    speciesMFP[i] = 1.0/speciesMFP[i];

                    speciesMCT[i] = 1.0/speciesMCT[i];

                    //Bird, eq (4.77)
                    MFP += speciesMFP[i]*nDensP/totalNumberDensity;

                    //Bird, eq (1.38)
                    MCT += speciesMCT[i]*nDensP/totalNumberDensity;
                }
            }

            // Calculate time-step to mean collision time ratio
            const scalar deltaT = mesh_.time().deltaTValue();

            scalar timeStepMCTRatio = deltaT/MCT;

            // Calculate cell size to mean free path ratio
            scalar largestCellDimension = 0.0;

            point minPoint = vector(GREAT, GREAT, GREAT);
            point maxPoint = vector(-GREAT, -GREAT, -GREAT);
            const List<label>& cellNodes = mesh_.cellPoints()[cell];

            forAll(cellNodes, node) 
            {
                const point& cellPoint = mesh_.points()[cellNodes[node]];
                minPoint.x() = min(minPoint.x(),cellPoint.x());
                minPoint.y() = min(minPoint.y(),cellPoint.y());
                minPoint.z() = min(minPoint.z(),cellPoint.z());
                maxPoint.x() = max(maxPoint.x(),cellPoint.x());
                maxPoint.y() = max(maxPoint.y(),cellPoint.y());
                maxPoint.z() = max(maxPoint.z(),cellPoint.z());                
            }

            vector cellSizeMFPRatio = (maxPoint-minPoint)/MFP;

            if (cloud_.cellCollModel(cell) == cloud_.binCollModel())
            {
                
                forAll(solutionDimensions, dim)
                {
                    if (solutionDimensions[dim])
                    {
                        cloud_.subcellLevels()[cell][dim] = 
                            label(min(cloud_.dynamicAdapter().maxSubcellLevels(),
                            max(cloud_.dynamicAdapter().minSubcellLevels(),cellSizeMFPRatio[dim]/cloud_.dynamicAdapter().maxSubcellSizeMFPRatio()))+0.5);
                    }
                    else
                    {
                        cloud_.subcellLevels()[cell][dim] = label(1.0);
                    }
                }
            }
            else
            {
                forAll(solutionDimensions, dim)
                {
                    if (solutionDimensions[dim])
                    {
                        cloud_.subcellLevels()[cell][dim] = cloud_.dynamicAdapter().minSubcellLevels();
                    }
                    else
                    {
                        cloud_.subcellLevels()[cell][dim] = label(1.0);
                    }
                }   
                 
            }
        }
        cloud_.subcellLevels().correctBoundaryConditions();
    }
    

    //Compute cell weights in case of cell weighted simulation    
    if (cloud_.cellWeighted())
    {
        forAll(mesh_.cells(), cell)
        {

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
