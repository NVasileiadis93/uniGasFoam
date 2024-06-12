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

#include "uniGasDynamicAdapter.H"
#include "uniGasCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
uniGasDynamicAdapter::uniGasDynamicAdapter
(
    const dictionary& dict,
    Time& time,
    fvMesh& mesh,
    uniGasCloud& cloud
)
:
    dict_(dict),
    time_(time),
    mesh_(mesh),
    cloud_(cloud),
    minSubcellLevels_(1),
    maxSubcellLevels_(10),
    timeStepAdaptation_(false),
    subcellAdaptation_(false),
    cellWeightAdaptation_(false),
    adaptationInterval_(),
    maxTimeStepMCTRatio_(),
    maxCourantNumber_(),
    maxSubcellSizeMFPRatio_(),
    smoothingPasses_(25),
    theta_(0.2),
    timeSteps_(0),
    timeAvCounter_(0),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    linearKEMeanXnParticle_(mesh_.nCells(), 0.0),
    momentumMeanXnParticle_(mesh_.nCells(), Zero),
    nParcelsXnParticle_(),
    rhoN_
    (
        IOobject
        (
            "rhoN",
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
            "translationalT",
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
            "UMean",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimVelocity, Zero),
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
    ),
    timeStepMCTRatio_
    (
        IOobject
        (
            "timeStepMCTRatio_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    courantNumber_
    (
        IOobject
        (
            "courantNumber_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ), 
    cellSizeMFPRatio_
    (
        IOobject
        (
            "cellSizeToMFPRatio",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    prevCellSizeMFPRatio_
    (
        IOobject
        (
            "prevCellSizeMFPRatio",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    )
{

    typeIds_.setSize(cloud_.typeIdList().size());
    forAll(typeIds_,iD)
    {
        typeIds_[iD] = iD;
    }

    nParcelsXnParticle_.setSize(typeIds_.size());
    for (auto& f : nParcelsXnParticle_)
    {
        f.setSize(mesh_.nCells());
    }

    if (cloud_.adaptive())
    {

        dictionary adaptationDict = dict.subDict("adaptiveProperties");

        timeStepAdaptation_ = adaptationDict.getOrDefault<bool>("timeStepAdaptation",false);
        
        subcellAdaptation_ = adaptationDict.getOrDefault<bool>("subcellAdaptation",false);

        cellWeightAdaptation_ = adaptationDict.getOrDefault<bool>("cellWeightAdaptation",false);
        
        adaptationInterval_ = adaptationDict.get<label>("adaptationInterval");

        if (timeStepAdaptation_)
        {
            maxTimeStepMCTRatio_ = adaptationDict.getOrDefault<scalar>("maxTimeStepMCTRatio", 0.2);
            maxCourantNumber_ = adaptationDict.getOrDefault<scalar>("maxCourantNumber", 0.5);
        }
        if (subcellAdaptation_)
        {
            maxSubcellSizeMFPRatio_ = adaptationDict.getOrDefault<scalar>("maxSubcellSizeMFPRatio",0.5);
        }

    }

}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uniGasDynamicAdapter::calculateAdaptationQuantities
(
    const label& cell,
    const scalar& rhoN,
    const scalar& transT,
    const vector& U,
    const scalarList& speciesRhoN,
    scalar& timeStepMCTRatio,
    vector& courantNumber,
    vector& cellSizeMFPRatio
)
{

    if (transT > SMALL)
    {

        scalar u0 = 0.0;
        scalarList speciesMCR(typeIds_.size(), 0.0);
        scalarList speciesMFP(typeIds_.size(), 0.0);

        forAll(typeIds_, iD)
        {
            
            for (label qspec=0; qspec<typeIds_.size(); ++qspec)
            {
                scalar dPQ = 0.5*(cloud_.constProps(typeIds_[iD]).d() + cloud_.constProps(typeIds_[qspec]).d());

                scalar omegaPQ = 0.5*(cloud_.constProps(typeIds_[iD]).omega() + cloud_.constProps(typeIds_[qspec]).omega());

                scalar massRatio = cloud_.constProps(typeIds_[iD]).mass()/cloud_.constProps(typeIds_[qspec]).mass();

                if (speciesRhoN[qspec] > VSMALL && transT > VSMALL)
                {

                    scalar reducedMass = 
                        cloud_.constProps(typeIds_[iD]).mass()*cloud_.constProps(typeIds_[qspec]).mass()
                        /(cloud_.constProps(typeIds_[iD]).mass() + cloud_.constProps(typeIds_[qspec]).mass());

                    // //Bird, eq (4.74)
                    speciesMCR[iD] += 
                        (2.0*sqrt(pi)*dPQ*dPQ*speciesRhoN[qspec]*pow(transT/cloud_.collTref(),1.0 - omegaPQ)*sqrt(2.0*physicoChemical::k.value()*cloud_.collTref()/reducedMass)); 

                    //Bird, eq (4.76)
                    speciesMFP[iD] += pi*dPQ*dPQ*speciesRhoN[qspec]*pow(cloud_.collTref()/transT,omegaPQ - 0.5)*sqrt(1.0 + massRatio);

                }
            }

            if (speciesRhoN[iD] > VSMALL)
            {
                u0 = max(u0,std::sqrt(2.0*Foam::constant::physicoChemical::k.value()/cloud_.constProps(typeIds_[iD]).mass()*transT));
            }

        }

        scalar MCR = 0.0;
        scalar MFP = 0.0;
        forAll(speciesMFP, iD)
        {
            if (speciesRhoN[iD] > VSMALL)
            {

                speciesMFP[iD] = 1.0/speciesMFP[iD];

                //Bird, eq (1.38)
                MCR += speciesMCR[iD]*speciesRhoN[iD]/rhoN;

                //Bird, eq (4.77)
                MFP += speciesMFP[iD]*speciesRhoN[iD]/rhoN;

            }
        }

        // Calculate cell dimensions
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

        // Calculate time-step to mean collision time ratio
        scalar deltaT = mesh_.time().deltaTValue();
        timeStepMCTRatio = deltaT*MCR;

        // Calculate Courant number
        courantNumber.x() = max(u0,U.x())*deltaT/(maxPoint.x()-minPoint.x());
        courantNumber.y() = max(u0,U.y())*deltaT/(maxPoint.y()-minPoint.y());
        courantNumber.z() = max(u0,U.z())*deltaT/(maxPoint.z()-minPoint.z());

        // Calculate cell size to mean free path ratio
        cellSizeMFPRatio = (maxPoint-minPoint)/MFP;

    }
    else
    {
        timeStepMCTRatio = 0.0;
        courantNumber = vector::zero;
        cellSizeMFPRatio = prevCellSizeMFPRatio_[cell];
    }
    
}

void uniGasDynamicAdapter::calculateTimeStep()
{

    // Calculate maximum time step to mean collision time ratio
    scalar instMaxTimeStepMCTRatio = 0e0;
    forAll(mesh_.cells(), cell)
    {
        if (cloud_.cellCollModel(cell) == cloud_.binCollModel() && instMaxTimeStepMCTRatio < timeStepMCTRatio_[cell])
        { 
            instMaxTimeStepMCTRatio = timeStepMCTRatio_[cell];
        }
    }

    // Calculate maximum Courant number
    scalar instMaxCourantNumber = 0e0;
    forAll(mesh_.cells(), cell)
    {
        forAll(cloud_.solutionDimensions(), dim)
        {
            if (cloud_.solutionDimensions()[dim] && instMaxCourantNumber < courantNumber_[cell][dim])
            {
                instMaxCourantNumber = courantNumber_[cell][dim];
            }
        }
    }

    // Set time-step
    scalar deltaT = mesh_.time().deltaTValue();
    if (instMaxTimeStepMCTRatio > VSMALL && instMaxCourantNumber < VSMALL)
    {
        deltaT *= maxTimeStepMCTRatio_/instMaxTimeStepMCTRatio;
    }
    else if (instMaxTimeStepMCTRatio < VSMALL && instMaxCourantNumber > VSMALL)
    {
        deltaT *= maxCourantNumber_/instMaxCourantNumber;
    }
    else if (instMaxTimeStepMCTRatio > VSMALL && instMaxCourantNumber > VSMALL)
    {
        deltaT *= min(maxTimeStepMCTRatio_/instMaxTimeStepMCTRatio,maxCourantNumber_/instMaxCourantNumber);
    }

    // Find minimum time step  across all processors
    if (Pstream::parRun())
    {
        reduce(deltaT, minOp<scalar>());
    }
    time_.setDeltaT(deltaT);

    //Delete for release
    label nUniGasParticles = cloud_.size();
    reduce(nUniGasParticles, sumOp<label>());
    if (Pstream::myProcNo() == 0 )
    {
        std::ofstream outfile;
        outfile.open("timeStepHistory.log", std::ios_base::app);
        outfile << mesh_.time().timeName() << " " << nUniGasParticles << " " << deltaT << "\n";
    }

}  

vector uniGasDynamicAdapter::calculateSubcellLevels
(
    const label& cell,
    const vector& cellSizeMFPRatio
)
{

    vector subcellLevels;

    if (cloud_.cellCollModel(cell) == cloud_.binCollModel())
    { 
        forAll(cloud_.solutionDimensions(), dim)
        {
            if (cloud_.solutionDimensions()[dim])
            {
                subcellLevels[dim] = label(min(maxSubcellLevels_,max(minSubcellLevels_,std::ceil(cellSizeMFPRatio[dim]/maxSubcellSizeMFPRatio_))));
            }
            else
            {
                subcellLevels[dim] = label(1.0);
            }
        }
    }
    else
    {
        forAll(cloud_.solutionDimensions(), dim)
        {
            if (cloud_.solutionDimensions()[dim])
            {
                subcellLevels[dim] = label(minSubcellLevels_);
            }
            else
            {
                subcellLevels[dim] = label(1.0);
            }
        }    
    }  

    return subcellLevels;

}

scalar uniGasDynamicAdapter::calculateCellWeightFactor
(
    const label& cell,
    const scalar& rhoN
)
{

    scalar RWF = cloud_.axiRWF(mesh_.C()[cell]);
    const vector& subcellLevels = cloud_.subcellLevels()[cell];
    const scalar nSubcells = subcellLevels.x()*subcellLevels.y()*subcellLevels.z();

    return (rhoN*mesh_.V()[cell])/(cloud_.particlesPerSubcell()*nSubcells*cloud_.nParticle()*RWF);

}  

void uniGasDynamicAdapter::smoothCellWeightFactor
(
    volScalarField& rhoN,
    volScalarField& cellWeightFactor
)
{

    //Initial smoothing
    fvc::smooth(cellWeightFactor, 1.3);

    scalar maxCellWeightRatio;
    label smoothingPasses = 0;
    do
    {

        smoothingPasses++;

        cellWeightFactor = fvc::average(fvc::interpolate(cellWeightFactor));
        cellWeightFactor.correctBoundaryConditions(); 

        forAll(mesh_.cells(), cell)
        {

            scalar RWF = cloud_.axiRWF(mesh_.C()[cell]);
            const vector& subcellLevels = cloud_.subcellLevels()[cell];
            const scalar nSubcells = subcellLevels.x()*subcellLevels.y()*subcellLevels.z();
            cellWeightFactor[cell] = 
                max(min((rhoN[cell]*mesh_.V()[cell])/(cloud_.minParticlesPerSubcell()*nSubcells*cloud_.nParticle()*RWF),cellWeightFactor[cell]),SMALL);

        }
        cellWeightFactor.correctBoundaryConditions();

        maxCellWeightRatio = VSMALL;
        forAll(mesh_.faces(), face)
        {

            if (mesh_.isInternalFace(face))
            {

                scalar ownerCWF = cellWeightFactor[mesh_.faceOwner()[face]];
                scalar neighbourCWF = cellWeightFactor[mesh_.faceNeighbour()[face]];
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

void uniGasDynamicAdapter::adapt()
{

    const scalar& deltaT = mesh_.time().deltaTValue();

    timeSteps_++;

    timeAvCounter_ += deltaT;

    // get cell measurements
    auto& cm = cloud_.cellPropMeasurements();

    forAll(cm.rhoNMean(), iD)
    {

        rhoNMean_ += deltaT*cm.rhoNMean()[iD];

        rhoNMeanXnParticle_ += deltaT*cm.rhoNMeanXnParticle()[iD];
        rhoMMeanXnParticle_ += deltaT*cm.rhoMMeanXnParticle()[iD];
        momentumMeanXnParticle_ += deltaT*cm.momentumMeanXnParticle()[iD];
        linearKEMeanXnParticle_ += deltaT*cm.linearKEMeanXnParticle()[iD];

        nParcelsXnParticle_[iD] += deltaT*cm.nParcelsXnParticle()[iD];

    }

    if (timeSteps_ == adaptationInterval_)
    {

        const auto& meshV = mesh_.V();

        // Computing internal fields
        scalar minRhoN = VGREAT;
        scalar minTranslationalT = VGREAT;
        forAll(rhoNMean_, cell)
        {

            if (rhoNMean_[cell] > VSMALL)
            {

                rhoN_[cell] = rhoNMeanXnParticle_[cell]/(meshV[cell]*timeAvCounter_);

                scalar rhoMMean = rhoMMeanXnParticle_[cell]/(meshV[cell]*timeAvCounter_);
                UMean_[cell] = momentumMeanXnParticle_[cell]/(rhoMMean*meshV[cell]*timeAvCounter_);

                scalar linearKEMean = 0.5*linearKEMeanXnParticle_[cell]/(meshV[cell]*timeAvCounter_);
                scalar rhoNMean = rhoNMeanXnParticle_[cell]/(meshV[cell]*timeAvCounter_);
                translationalT_[cell] =
                    2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                   *(
                        linearKEMean
                      - 0.5*rhoMMean*(UMean_[cell] & UMean_[cell])
                    );

                if (translationalT_[cell] > VSMALL)
                {
                    if (minRhoN >  rhoN_[cell])
                    {
                        minRhoN = rhoN_[cell];
                    }
                    if (minTranslationalT >  translationalT_[cell])
                    {
                        minTranslationalT = translationalT_[cell];
                    }
                }
            }
            else
            {
                rhoN_[cell] = 0.0;
                translationalT_[cell] = 0.0;
                UMean_[cell] = vector::zero;
            }

        }


        // Calculate cell size to mean free path ratio
        forAll(mesh_.cells(), cell)
        {

            scalarList speciesRhoN(typeIds_.size(), 0.0);
            forAll(typeIds_,iD)
            {
                speciesRhoN[iD] = nParcelsXnParticle_[iD][cell]/(meshV[cell]*timeAvCounter_);
            }

            // set zero density and temperature to the smallest non-zero ones to avoid errors in cell weighting factors
            if (rhoN_[cell] == 0.0)
            {
                rhoN_[cell] = minRhoN;
                forAll(typeIds_,iD)
                {
                    speciesRhoN[iD] = rhoN_[cell]/typeIds_.size();
                }
            }
            if (translationalT_[cell] == 0.0)
            {
                translationalT_[cell] = minTranslationalT;
            }          

            calculateAdaptationQuantities
                (
                    cell,
                    rhoN_[cell], 
                    translationalT_[cell],
                    UMean_[cell],
                    speciesRhoN,
                    timeStepMCTRatio_[cell],
                    courantNumber_[cell],
                    cellSizeMFPRatio_[cell]
                );    
        
        }
        timeStepMCTRatio_.correctBoundaryConditions();
        courantNumber_.correctBoundaryConditions();
        cellSizeMFPRatio_.correctBoundaryConditions();

        // Smooth macroscopic fields
        for (label pass = 1; pass <= smoothingPasses_; ++pass)
        {
            cellSizeMFPRatio_ = fvc::average(fvc::interpolate(cellSizeMFPRatio_));
            cellSizeMFPRatio_.correctBoundaryConditions();
        }

        // Adapt time step
        if (timeStepAdaptation_)
        {                      
            calculateTimeStep();
        }

        // Adapt subcell levels
        if (subcellAdaptation_)
        {
            forAll(mesh_.cells(), cell)
            {
                cloud_.subcellLevels()[cell] = calculateSubcellLevels
                                                (
                                                    cell,
                                                    cellSizeMFPRatio_[cell]
                                                );
            }                           
            cloud_.subcellLevels().correctBoundaryConditions();
        }

        // Adapt cell weight factor
        if (cellWeightAdaptation_)
        {

            // Calculate cell weight factor
            forAll(mesh_.cells(), cell)
            {
                cellWeightFactor_[cell] = calculateCellWeightFactor
                                            (
                                                cell,
                                                rhoN_[cell]
                                            );
            }                           
            cloud_.subcellLevels().correctBoundaryConditions();

            // Smooth cell weight factor
            smoothCellWeightFactor
                (
                    rhoN_,
                    cellWeightFactor_
                );

            // Time average cell weight factor
            cloud_.cellWeightFactor() = theta_*cellWeightFactor_ + (1.0-theta_)*cloud_.cellWeightFactor();
            cellWeightFactor_ = cloud_.cellWeightFactor();

        }

        // Reset
        timeSteps_ = 0;
        
        timeAvCounter_ = 0.0;

        forAll(rhoN_, cell)
        {

            rhoNMean_[cell] = 0.0;
            rhoNMeanXnParticle_[cell] = 0.0;
            rhoMMeanXnParticle_[cell] = 0.0;
            momentumMeanXnParticle_[cell] = vector::zero;
            linearKEMeanXnParticle_[cell] = 0.0;

            forAll(typeIds_, iD)
            {
                nParcelsXnParticle_[iD][cell] = 0.0;
            }
        }

        // Store old data
        prevCellSizeMFPRatio_ = cellSizeMFPRatio_;

    }

}

void uniGasDynamicAdapter::setInitialConfiguration
(
    const Field<scalar>& speciesRhoN,
    const scalar& transT,
    const vector& U
)
{

    // Calculate total number density
    scalar rhoN = 0.0;
    forAll(typeIds_,iD)
    {
        rhoN += speciesRhoN[iD];
    }

    // Calculate cell size to mean free path ratio
    forAll(mesh_.cells(), cell)
    {
        calculateAdaptationQuantities
            (
                cell,
                rhoN, 
                transT, 
                U,
                speciesRhoN,
                timeStepMCTRatio_[cell],
                courantNumber_[cell],
                cellSizeMFPRatio_[cell]
            );
    }
    timeStepMCTRatio_.correctBoundaryConditions();
    courantNumber_.correctBoundaryConditions();
    cellSizeMFPRatio_.correctBoundaryConditions();

    // Smooth cell size to mean free path ratio
    for (label pass = 1; pass <= smoothingPasses_; ++pass)
    {
        timeStepMCTRatio_ = fvc::average(fvc::interpolate(timeStepMCTRatio_));
        courantNumber_ = fvc::average(fvc::interpolate(courantNumber_));
        cellSizeMFPRatio_ = fvc::average(fvc::interpolate(cellSizeMFPRatio_));
        timeStepMCTRatio_.correctBoundaryConditions();
        courantNumber_.correctBoundaryConditions();
        cellSizeMFPRatio_.correctBoundaryConditions();
    }

    // Store old data
    prevCellSizeMFPRatio_ = cellSizeMFPRatio_;

    // Adapt time step
    if (timeStepAdaptation_)
    {                      
        calculateTimeStep();
    }

    // Adapt subcell levels
    if (subcellAdaptation_)
    {
        forAll(mesh_.cells(), cell)
        {
            cloud_.subcellLevels()[cell] = calculateSubcellLevels
                                            (
                                                cell,
                                                cellSizeMFPRatio_[cell]
                                            );
        }                           
        cloud_.subcellLevels().correctBoundaryConditions();
    }

}

void uniGasDynamicAdapter::setInitialConfiguration
(
    const cellZone& zone,
    const Field<scalar>& speciesRhoN,
    const scalar& transT,
    const vector& U
)
{

    // Calculate total number density
    scalar rhoN = 0.0;
    forAll(typeIds_,iD)
    {
        rhoN += speciesRhoN[iD];
    }

    // Calculate cell size to mean free path ratio
    if (zone.size())
    {
        for (const label cell : zone)
        {

            calculateAdaptationQuantities
              (
                  cell,
                  rhoN, 
                  transT,
                  U,
                  speciesRhoN,
                  timeStepMCTRatio_[cell],
                  courantNumber_[cell],
                  cellSizeMFPRatio_[cell]
              );
        }
    }
    timeStepMCTRatio_.correctBoundaryConditions();
    courantNumber_.correctBoundaryConditions();
    cellSizeMFPRatio_.correctBoundaryConditions();

    // Store old data
    prevCellSizeMFPRatio_ = cellSizeMFPRatio_;

    // Adapt time step
    if (timeStepAdaptation_)
    {                      
        calculateTimeStep();
    }

    // Adapt subcell levels
    if (subcellAdaptation_)
    {
        if (zone.size())
        {
            for (const label cell : zone)
            {
                cloud_.subcellLevels()[cell] = calculateSubcellLevels
                                                (
                                                    cell,
                                                    cellSizeMFPRatio_[cell]
                                                );
            }            
        }               
        cloud_.subcellLevels().correctBoundaryConditions();
    }

}

void uniGasDynamicAdapter::setInitialConfiguration
(
    const List<autoPtr<volScalarField>>& speciesRhoNPtr,
    const volScalarField& transT,
    const volVectorField& U
)
{

    // Calculate total number density
    volScalarField rhoN
    (
        IOobject
        (
            "rhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, 0.0)
    );

    forAll(typeIds_,iD)
    {
        const volScalarField& speciesRhoN = speciesRhoNPtr[iD];
        rhoN += speciesRhoN;
    }

    // Calculate cell size to mean free path ratio
    forAll(mesh_.cells(), cell)
    {

        scalarList speciesRhoN(typeIds_.size(), 0.0);
        forAll(typeIds_,iD)
        {
            const volScalarField& speciesRhoNID = speciesRhoNPtr[iD];
            speciesRhoN[iD] = speciesRhoNID[cell];
        }

        calculateAdaptationQuantities
        (
            cell,
            rhoN[cell], 
            transT[cell],
            U[cell],
            speciesRhoN,
            timeStepMCTRatio_[cell],
            courantNumber_[cell],
            cellSizeMFPRatio_[cell]
        );


    }
    timeStepMCTRatio_.correctBoundaryConditions();
    courantNumber_.correctBoundaryConditions();
    cellSizeMFPRatio_.correctBoundaryConditions();

    // Store old data
    prevCellSizeMFPRatio_ = cellSizeMFPRatio_;

    // Adapt time step
    if (timeStepAdaptation_)
    {                      
        calculateTimeStep();
    }

    // Adapt subcell levels
    if (subcellAdaptation_)
    {
        forAll(mesh_.cells(), cell)
        {
            cloud_.subcellLevels()[cell] = calculateSubcellLevels
                                            (
                                                cell,
                                                cellSizeMFPRatio_[cell]
                                            );
        }                           
        cloud_.subcellLevels().correctBoundaryConditions();
    }

}

} // End namespace Foam

// ************************************************************************* //

