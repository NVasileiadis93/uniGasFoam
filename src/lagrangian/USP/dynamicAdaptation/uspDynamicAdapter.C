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
    minSubcellLevels_(2),
    maxSubcellLevels_(5),
    cellWeightAdaptation_(false),
    subcellAdaptation_(false),
    adaptationInterval_(),
    maxSubcellSizeMFPRatio_(),
    smoothingPasses_(10),
    theta_(0.1),
    timeSteps_(0),
    nAvTimeSteps_(0),
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
    cellSizeMFPRatio_
    (
        IOobject
        (
            "cellSizeToMFPRatioAdapt",
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
            "prevCellSizeMFPRatio_",
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
        subcellAdaptation_ = adaptationDict.getOrDefault<bool>("subcellAdaptation",false);
        cellWeightAdaptation_ = adaptationDict.getOrDefault<bool>("cellWeightAdaptation",false);
        if (subcellAdaptation_)
        {
            maxSubcellSizeMFPRatio_ = adaptationDict.get<scalar>("maxSubcellSizeMFPRatio");
        }
        adaptationInterval_ = adaptationDict.get<label>("adaptationInterval");
    }

}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector uspDynamicAdapter::calculateCellSizeMFPRatio
(
    const label& cell,
    const scalar& rhoN,
    const scalar& transT,
    const scalarList& speciesRhoN
)
{

    if (transT > SMALL)
    {

        scalarList speciesMFP(typeIds_.size(), 0.0);

        forAll(typeIds_, iD)
        {
            label qspec = 0;

            for (qspec=0; qspec<typeIds_.size(); ++qspec)
            {
                scalar dPQ = 0.5*(cloud_.constProps(typeIds_[iD]).d() + cloud_.constProps(typeIds_[qspec]).d());

                scalar omegaPQ = 0.5*(cloud_.constProps(typeIds_[iD]).omega() + cloud_.constProps(typeIds_[qspec]).omega());

                scalar massRatio = cloud_.constProps(typeIds_[iD]).mass()/cloud_.constProps(typeIds_[qspec]).mass();

                if (speciesRhoN[qspec] > VSMALL && transT > VSMALL)
                {

                    //Bird, eq (4.76)
                    speciesMFP[iD] += pi*dPQ*dPQ*speciesRhoN[qspec]*pow(cloud_.collTref()/transT,omegaPQ - 0.5)*sqrt(1.0 + massRatio);

                }
            }

        }

        scalar MFP = 0.0;
        forAll(speciesMFP, iD)
        {
            if (speciesRhoN[iD] > VSMALL)
            {

                speciesMFP[iD] = 1.0/speciesMFP[iD];

                //Bird, eq (4.77)
                MFP += speciesMFP[iD]*speciesRhoN[iD]/rhoN;

            }
        }

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

        return (maxPoint-minPoint)/MFP;

    }
    else
    {
        return prevCellSizeMFPRatio_[cell];
    }
    
}

vector uspDynamicAdapter::calculateSubcellLevels
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
                subcellLevels[dim] = label(min(maxSubcellLevels_,max(minSubcellLevels_,cellSizeMFPRatio[dim]/maxSubcellSizeMFPRatio_))+0.5);
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
                subcellLevels[dim] = minSubcellLevels_;
            }
            else
            {
                subcellLevels[dim] = label(1.0);
            }
        }    
    }

    return subcellLevels;

}

scalar uspDynamicAdapter::calculateCellWeightFactor
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

void uspDynamicAdapter::smoothCellWeightFactor
(
    volScalarField& rhoN,
    volScalarField& cellWeightFactor
)
{

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

void uspDynamicAdapter::adapt()
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

        nParcelsXnParticle_[iD] += cm.nParcelsXnParticle()[iD];

    }

    if (timeSteps_ == adaptationInterval_)
    {

        const auto& meshCC = mesh_.cellCentres();
        const auto& meshV = mesh_.V();

        // Computing internal fields
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
                speciesRhoN[iD] = nParcelsXnParticle_[iD][cell]/(meshV[cell]*nAvTimeSteps_);
            }

            cellSizeMFPRatio_[cell] = calculateCellSizeMFPRatio
                                        (
                                            cell,
                                            rhoN_[cell], 
                                            translationalT_[cell], 
                                            speciesRhoN
                                        );
        }
        cellSizeMFPRatio_.correctBoundaryConditions();

        // Smooth cell size to mean free path ratio
        //for (label pass = 1; pass <= smoothingPasses_; ++pass)
        //{
        //    cellSizeMFPRatio_ = fvc::average(fvc::interpolate(cellSizeMFPRatio_));
        //    cellSizeMFPRatio_.correctBoundaryConditions();
        //}

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
        
        nAvTimeSteps_ = 0;

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

void uspDynamicAdapter::setInitialConfiguration
(
    const Field<scalar>& speciesRhoN,
    const scalar& transT
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
            cellSizeMFPRatio_[cell] = calculateCellSizeMFPRatio
                                        (
                                            cell,
                                            rhoN, 
                                            transT, 
                                            speciesRhoN
                                        );
        }
        cellSizeMFPRatio_.correctBoundaryConditions();

        // Store old size to mean free path ratio
        prevCellSizeMFPRatio_ = cellSizeMFPRatio_;

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

void uspDynamicAdapter::setInitialConfiguration
(
    const cellZone& zone,
    const Field<scalar>& speciesRhoN,
    const scalar& transT
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
            cellSizeMFPRatio_[cell] = calculateCellSizeMFPRatio
                                        (
                                            cell,
                                            rhoN, 
                                            transT, 
                                            speciesRhoN
                                        );
        }
    }
    cellSizeMFPRatio_.correctBoundaryConditions();

    // Store old size to mean free path ratio
    prevCellSizeMFPRatio_ = cellSizeMFPRatio_;

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

void uspDynamicAdapter::setInitialConfiguration
(
    const List<autoPtr<volScalarField>>& speciesRhoNPtr,
    const volScalarField& transT
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

        cellSizeMFPRatio_[cell] = calculateCellSizeMFPRatio
                                    (
                                        cell,
                                        rhoN[cell], 
                                        transT[cell], 
                                        speciesRhoN
                                    );
    }
    cellSizeMFPRatio_.correctBoundaryConditions();

    // Store old size to mean free path ratio
    prevCellSizeMFPRatio_ = cellSizeMFPRatio_;

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

