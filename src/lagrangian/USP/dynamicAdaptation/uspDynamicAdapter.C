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

#include "uspDynamicAdapter.H"
#include "uspCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
uspDynamicAdapter::uspDynamicAdapter
(
    const dictionary& dict,
    const fvMesh& mesh,
    uspCloud& cloud
)
:
    dict_(dict),
    mesh_(mesh),
    cloud_(cloud),
    rndGen_(cloud.rndGen()),
    minSubcellLevels_(2),
    maxSubcellLevels_(5),
    timeStepAdaptation_(false),
    cellWeightAdaptation_(false),
    subcellAdaptation_(false),
    timeInterval_(),
    maxSubcellSizeMFPRatio_(),
    Tref_(),
    smoothingPasses_(10),
    theta_(0.8),
    timeSteps_(0),
    nAvTimeSteps_(0),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    linearKEMeanXnParticle_(mesh_.nCells(), 0.0),
    momentumMeanXnParticle_(mesh_.nCells(), Zero),
    MFP_(mesh_.nCells(), 0.0),
    MCT_(mesh_.nCells(), 0.0),
    nParcels_(),
    nParcelsXnParticle_(),
    rhoN_
    (
        IOobject
        (
            "rhoN_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    translationalT_
    (
        IOobject
        (
            "translationalT_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    UMean_
    (
        IOobject
        (
            "UMean_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimVelocity, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    cellSizeMFPRatio_
    (
        IOobject
        (
            "cellSizeToMFPRatioAdapt",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    timeStepMCTRatio_
    (
        IOobject
        (
            "timeStepToMCTRatio",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    cellWeightFactor_
    (
        IOobject
        (
            "cellWeightFactor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    )
{

    typeIds_.setSize(cloud_.typeIdList().size());
    forAll(typeIds_,iD)
    {
        typeIds_[iD] = iD;
    }

    nParcels_.setSize(typeIds_.size());

    for (auto& f : nParcels_)
    {
        f.setSize(mesh_.nCells());
    }

    nParcelsXnParticle_.setSize(typeIds_.size());
    for (auto& f : nParcelsXnParticle_)
    {
        f.setSize(mesh_.nCells());
    }

    mfp_.setSize(typeIds_.size());

    for (auto& f : mfp_)
    {
        f.setSize(mesh_.nCells());
    }

    mct_.setSize(typeIds_.size());

    for (auto& f : mct_)
    {
        f.setSize(mesh_.nCells());
    }

    if (cloud_.dynamicAdaptation())
    {
        timeStepAdaptation_ = dict.subDict("dynamicSimulationProperties").getOrDefault<bool>("timeStepAdaptation",false);
        subcellAdaptation_ = dict.subDict("dynamicSimulationProperties").getOrDefault<bool>("subcellAdaptation",false);
        cellWeightAdaptation_ = dict.subDict("dynamicSimulationProperties").getOrDefault<bool>("cellWeightAdaptation",false);
        if (subcellAdaptation_)
        {
            maxSubcellSizeMFPRatio_ = dict.subDict("dynamicSimulationProperties").get<scalar>("maxSubcellSizeMFPRatio");
        }
        timeInterval_ = dict.subDict("dynamicSimulationProperties").get<label>("timeInterval");
        Tref_ = dict.subDict("collisionProperties").get<scalar>("Tref");
    }

}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uspDynamicAdapter::calculateProperties()
{
    
    timeSteps_++;

    nAvTimeSteps_++;

    // get cell measurements
    auto& cm = cloud_.cellPropMeasurements();

    forAll(cm.rhoNMean(), iD)
    {

        rhoNMean_ += cm.rhoNMean()[iD];

        rhoNMeanXnParticle_ += cm.rhoNMeanXnParticle()[iD];
        rhoMMeanXnParticle_ += cm.rhoMMeanXnParticle()[iD];
        momentumMeanXnParticle_ += cm.momentumMeanXnParticle()[iD];
        linearKEMeanXnParticle_ += cm.linearKEMeanXnParticle()[iD];

        nParcels_[iD] += cm.nParcels()[iD];
        nParcelsXnParticle_[iD] += cm.nParcelsXnParticle()[iD];

    }

    if (timeSteps_ == timeInterval_)
    {

        const auto& meshCC = cloud_.mesh().cellCentres();
        const auto& meshV = cloud_.mesh().V();

        // computing internal fields
        forAll(rhoNMean_, cell)
        {
            if (rhoNMean_[cell] > VSMALL)
            {

                rhoN_[cell] = rhoNMeanXnParticle_[cell]/(meshV[cell]*nAvTimeSteps_);

                scalar rhoMMean = rhoMMeanXnParticle_[cell]/(meshV[cell]*nAvTimeSteps_);
                UMean_[cell] = momentumMeanXnParticle_[cell]/(rhoMMean*meshV[cell]*nAvTimeSteps_);

                scalar linearKEMean = 0.5*linearKEMeanXnParticle_[cell]/(meshV[cell]*nAvTimeSteps_);
                scalar rhoNMean = rhoNMeanXnParticle_[cell]/(meshV[cell]*nAvTimeSteps_);
                translationalT_[cell] =
                    2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                   *(
                        linearKEMean
                      - 0.5*rhoMMean*(UMean_[cell] & UMean_[cell])
                    );

                forAll(typeIds_, iD)
                {
                    label qspec = 0;

                    for (qspec=0; qspec<typeIds_.size(); ++qspec)
                    {
                        scalar dPQ = 0.5*(cloud_.constProps(typeIds_[iD]).d() + cloud_.constProps(typeIds_[qspec]).d());

                        scalar omegaPQ = 0.5*(cloud_.constProps(typeIds_[iD]).omega() + cloud_.constProps(typeIds_[qspec]).omega());

                        scalar massRatio = cloud_.constProps(typeIds_[iD]).mass()/cloud_.constProps(typeIds_[qspec]).mass();

                        if (nParcels_[qspec][cell] > VSMALL && translationalT_[cell] > VSMALL)
                        {
                            scalar nDensQ = (nParcelsXnParticle_[qspec][cell])/(meshV[cell]*nAvTimeSteps_);

                            scalar reducedMass = 
                                cloud_.constProps(typeIds_[iD]).mass()*cloud_.constProps(typeIds_[qspec]).mass()
                                /(cloud_.constProps(typeIds_[iD]).mass()+ cloud_.constProps(typeIds_[qspec]).mass());

                            //Bird, eq (4.76)
                            mfp_[iD][cell] += pi*dPQ*dPQ*nDensQ*pow(Tref_/translationalT_[cell],omegaPQ - 0.5)*sqrt(1.0 + massRatio);

                            // //Bird, eq (4.74)
                            mct_[iD][cell] +=
                                2.0*sqrt(pi)*dPQ*dPQ*nDensQ*pow(translationalT_[cell]/Tref_,1.0 - omegaPQ)*sqrt(2.0*physicoChemical::k.value()*Tref_/reducedMass); 
                        }
                    }

                }

                MFP_ = 0.0;
                MCT_ = 0.0;
                forAll(mfp_, iD)
                {
                    if (rhoN_[cell] > VSMALL)
                    {
                        scalar nDensP = (nParcelsXnParticle_[iD][cell])/(meshV[cell]*nAvTimeSteps_);

                        mfp_[iD][cell] = 1.0/mfp_[iD][cell];

                        mct_[iD][cell] = 1.0/mct_[iD][cell];

                        //Bird, eq (4.77)
                        MFP_[cell] += mfp_[iD][cell]*nDensP/rhoN_[cell];

                        //Bird, eq (1.38)
                        MCT_[cell] += mct_[iD][cell]*nDensP/rhoN_[cell];
                    }
                }

                // Calculate time-step to mean collision time ratio
                const scalar deltaT = mesh_.time().deltaTValue();

                timeStepMCTRatio_[cell] = deltaT/MCT_[cell];

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

                cellSizeMFPRatio_[cell] = (maxPoint-minPoint)/MFP_[cell];

            }
        }
        cellSizeMFPRatio_.correctBoundaryConditions();

        // smooth fields
        for (label pass=0; pass<smoothingPasses_; ++pass)
        {
            cellSizeMFPRatio_ = fvc::average(fvc::interpolate(cellSizeMFPRatio_));
            cellSizeMFPRatio_.correctBoundaryConditions();
        }
    }

}

void uspDynamicAdapter::update()
{

    if (timeSteps_ == timeInterval_)
    {

        const auto& meshCC = cloud_.mesh().cellCentres();
        const auto& meshV = cloud_.mesh().V();

        // Adapt time-step 
        if (timeStepAdaptation_)
        {
            //not coded yet
        }

        // Adapt subcell levels
        if (subcellAdaptation_)
        {
            const boolVector& solutionDimensions = cloud_.solutionDimensions(); 

            forAll(mesh_.cells(), cell)
            {
                if (cloud_.cellCollModel(cell) == cloud_.binCollModel())
                {
                    
                    forAll(solutionDimensions, dim)
                    {
                        if (solutionDimensions[dim])
                        {
                            cloud_.subcellLevels()[cell][dim] = label(min(maxSubcellLevels_,max(minSubcellLevels_,cellSizeMFPRatio_[cell][dim]/maxSubcellSizeMFPRatio_))+0.5);
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
                            cloud_.subcellLevels()[cell][dim] = minSubcellLevels_;
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

        // update cell weighting factor        
        if (cellWeightAdaptation_)
        {
            forAll(mesh_.cells(), cell)
            {
                scalar RWF = cloud_.axiRWF(meshCC[cell]);
                const vector& subcellLevels = cloud_.subcellLevels()[cell];
                const scalar nSubcells = subcellLevels.x()*subcellLevels.y()*subcellLevels.z();
                cellWeightFactor_[cell] = (rhoN_[cell]*meshV[cell])/(cloud_.particlesPerSubcell()*nSubcells*cloud_.nParticle()*RWF);
            }
            cellWeightFactor_.correctBoundaryConditions();

            // smooth cell weighting factor
            scalar maxCellWeightRatio;
            label smoothingPasses = 0;
            do
            {

                smoothingPasses++;

                cellWeightFactor_ = fvc::average(fvc::interpolate(cellWeightFactor_));
                cellWeightFactor_.correctBoundaryConditions(); 

                forAll(mesh_.cells(), cell)
                {

                    scalar RWF = cloud_.axiRWF(meshCC[cell]);
                    const vector& subcellLevels = cloud_.subcellLevels()[cell];
                    const scalar nSubcells = subcellLevels.x()*subcellLevels.y()*subcellLevels.z();
                    cellWeightFactor_[cell] = max((rhoN_[cell]*meshV[cell])/(cloud_.minParticlesPerSubcell()*nSubcells*cloud_.nParticle()*RWF),cellWeightFactor_[cell]);

                }
                cellWeightFactor_.correctBoundaryConditions();

                maxCellWeightRatio = VSMALL;
                forAll(mesh_.faces(), face)
                {

                    if (mesh_.isInternalFace(face))
                    {

                        scalar ownerCWF = cellWeightFactor_[mesh_.faceOwner()[face]];
                        scalar neighbourCWF = cellWeightFactor_[mesh_.faceNeighbour()[face]];
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

            // time average cell weight factor
            cloud_.cellWeightFactor() = theta_*cellWeightFactor_ + (1.0-theta_)*cloud_.cellWeightFactor();
            cellWeightFactor_ = cloud_.cellWeightFactor();

        }

    }

}

void uspDynamicAdapter::reset()
{

    if (timeSteps_ == timeInterval_)
    {

        // reset
        timeSteps_ = 0;

        forAll(rhoN_, cell)
        {

            nAvTimeSteps_ = 0;

            rhoNMean_[cell] = 0.0;
            rhoNMeanXnParticle_[cell] = 0.0;
            rhoMMeanXnParticle_[cell] = 0.0;
            momentumMeanXnParticle_[cell] = vector::zero;
            linearKEMeanXnParticle_[cell] = 0.0;

            forAll(typeIds_, iD)
            {
                nParcels_[iD][cell] = 0.0;
                nParcelsXnParticle_[iD][cell] = 0.0;
                mfp_[iD][cell] = 0.0;
                mct_[iD][cell] = 0.0;
            }
        }
    }

}

void uspDynamicAdapter::adapt()
{

    calculateProperties();

    update();

    reset();

}

}  // End namespace Foam

// ************************************************************************* //

