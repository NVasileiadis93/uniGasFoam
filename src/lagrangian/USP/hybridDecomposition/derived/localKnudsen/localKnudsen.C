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

#include "localKnudsen.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{

defineTypeNameAndDebug(localKnudsen, 0);

addToRunTimeSelectionTable(uspHybridDecomposition, localKnudsen, dictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::localKnudsen::localKnudsen
(
    const Time& t,
    const polyMesh& mesh,
    uspCloud& cloud
)
:
    uspHybridDecomposition(t, mesh, cloud),
    propsDict_(hybridDecompositionDict_.subDict(typeName + "Properties")),
    breakdownMax_(propsDict_.get<scalar>("breakdownMax")),
    theta_(propsDict_.getOrDefault<scalar>("theta",1.0)),
    smoothingPasses_(propsDict_.getOrDefault<scalar>("smoothingPasses",0)),    
    timeSteps_(0),
    timeAvCounter_(0.0),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    linearKEMeanXnParticle_(mesh_.nCells(), 0.0),
    momentumMeanXnParticle_(mesh_.nCells(), vector::zero),
    nParcelsXnParticle_(),
    KnGLL_
    (
        IOobject
        (
            "KnGLL",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    KnRho_
    (
        IOobject
        (
            "KnRho",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    KnT_
    (
        IOobject
        (
            "KnT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    KnU_
    (
        IOobject
        (
            "KnU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
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
    rhoM_
    (
        IOobject
        (
            "rhoM_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimVolume, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    p_
    (
        IOobject
        (
            "p_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure, Zero),
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
    )
{

    typeIds_.setSize(cloud_.typeIdList().size());
    forAll(typeIds_,iD)
    {
        typeIds_[iD] = iD;
    }

    nParcelsXnParticle_.setSize(typeIds_.size());

    for (auto& n : nParcelsXnParticle_)
    {
        n.setSize(mesh_.nCells(), 0.0);
    }

}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::localKnudsen::decompose()
{

    const scalar& deltaT = mesh_.time().deltaTValue();

    timeSteps_++;

    timeAvCounter_ += deltaT;

    // get cell measurements
    auto& cm = cloud_.cellPropMeasurements();

    forAll(cm.rhoNMean(), i)
    {

        rhoNMean_ += deltaT*cm.rhoNMean()[i];
        rhoMMean_  += deltaT*cm.rhoMMean()[i];
        linearKEMean_ += deltaT*cm.linearKEMean()[i];
        rhoNMeanXnParticle_ += deltaT*cm.rhoNMeanXnParticle()[i];
        rhoMMeanXnParticle_ += deltaT*cm.rhoMMeanXnParticle()[i];
        linearKEMeanXnParticle_ += deltaT*cm.linearKEMeanXnParticle()[i];
        momentumMeanXnParticle_ += deltaT*cm.momentumMeanXnParticle()[i];
        nParcelsXnParticle_[i] += deltaT*cm.nParcelsXnParticle()[i];

    }

    if (timeSteps_ == decompositionInterval_)
    {

        // computing internal fields
        forAll(rhoNMean_, cell)
        {
            if (rhoNMean_[cell] > VSMALL)
            {
                const scalar cellVolume = mesh_.cellVolumes()[cell];

                rhoN_[cell] = rhoNMeanXnParticle_[cell]/(timeAvCounter_*cellVolume);
                rhoM_[cell] = rhoMMeanXnParticle_[cell]/(timeAvCounter_*cellVolume);

                scalar rhoMMean = rhoMMeanXnParticle_[cell]/(cellVolume*timeAvCounter_);
                UMean_[cell] = momentumMeanXnParticle_[cell]/(rhoMMean*cellVolume*timeAvCounter_);

                scalar linearKEMean = 0.5*linearKEMeanXnParticle_[cell]/(cellVolume*timeAvCounter_);
                scalar rhoNMean = rhoNMeanXnParticle_[cell]/(cellVolume*timeAvCounter_);
                translationalT_[cell] =
                    2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                   *(
                        linearKEMean
                      - 0.5*rhoMMean*(UMean_[cell] & UMean_[cell])
                    );

                p_[cell] = rhoN_[cell]*physicoChemical::k.value()*translationalT_[cell];

            }
            else
            {
                rhoN_[cell] = 0.0;
                rhoM_[cell] = 0.0;
                p_[cell] = 0.0;
                translationalT_[cell] = 0.0;
                UMean_[cell] = vector::zero;
            }
        }

        rhoN_.correctBoundaryConditions();
        rhoM_.correctBoundaryConditions();
        p_.correctBoundaryConditions();
        translationalT_.correctBoundaryConditions();
        UMean_.correctBoundaryConditions();

        // smooth macroscopic quantities
        for (label pass=1; pass<=smoothingPasses_; pass++)
        {
            rhoM_ = fvc::average(fvc::interpolate(rhoM_));
            p_ = fvc::average(fvc::interpolate(p_));
            translationalT_ = fvc::average(fvc::interpolate(translationalT_));
            UMean_ = fvc::average(fvc::interpolate(UMean_)); 
            rhoM_.correctBoundaryConditions();
            p_.correctBoundaryConditions();
            translationalT_.correctBoundaryConditions();
            UMean_.correctBoundaryConditions();
        }

        scalarField maxMagGradRho(mesh_.nCells());
        scalarField maxMagGradT(mesh_.nCells());
        scalarField maxMagGradU(mesh_.nCells());

        forAll(mesh_.cells(), cell)
        {
        
            maxMagGradRho[cell] = 0.0;
            maxMagGradT[cell] = 0.0;
            maxMagGradU[cell] = 0.0;
            forAll(mesh_.cellCells()[cell], cellI)
            {
        
                label adjCell = mesh_.cellCells()[cell][cellI];
                scalar cellDistance = mag(mesh_.C()[adjCell]-mesh_.C()[cell]);
        
                maxMagGradRho[cell] = max(maxMagGradRho[cell], fabs(rhoM_[adjCell]-rhoM_[cell])/cellDistance);
                maxMagGradT[cell] = max(maxMagGradT[cell], fabs(translationalT_[adjCell]-translationalT_[cell])/cellDistance);
                maxMagGradU[cell] = max(maxMagGradU[cell], fabs(mag(UMean_[adjCell])-mag(UMean_[cell]))/cellDistance);
        
            }
        }        

        // Calculate KnGLL based on macroscopic quantities
        scalar instKnRho;
        scalar instKnT;
        scalar instKnU;
        scalar instKnGLL;

        forAll(mesh_.cells(), cell)
        {

            const scalar cellVolume = mesh_.cellVolumes()[cell];

            if (rhoNMean_[cell] > VSMALL && translationalT_[cell] > VSMALL)
            {

                scalarList speciesMFP(typeIds_.size(), 0.0);

                forAll(typeIds_, i)
                {
                    
                    for (label qspec=0; qspec<typeIds_.size(); ++qspec)
                    {
                        scalar dPQ = 0.5*(cloud_.constProps(i).d() + cloud_.constProps(qspec).d());

                        scalar omegaPQ = 0.5*(cloud_.constProps(i).omega() + cloud_.constProps(qspec).omega());

                        scalar massRatio = cloud_.constProps(i).mass()/cloud_.constProps(qspec).mass();

                        if (nParcelsXnParticle_[qspec][cell] > VSMALL)
                        {

                            scalar nDensQ = nParcelsXnParticle_[qspec][cell]/cellVolume;

                            //Bird, eq (4.76)
                            speciesMFP[i] += pi*dPQ*dPQ*nDensQ*pow(cloud_.collTref()/translationalT_[cell], omegaPQ-0.5)*sqrt(1.0 + massRatio); 

                        }
                    }

                }

                scalar MFP = 0.0;
                forAll(typeIds_, i)
                {
                    if (nParcelsXnParticle_[i][cell] > VSMALL)
                    {

                        speciesMFP[i] = 1.0/speciesMFP[i];

                        //Bird, eq (4.77)
                        MFP += speciesMFP[i]*nParcelsXnParticle_[i][cell]/(rhoN_[cell]*cellVolume);

                    }
                }

                scalar u0 = std::sqrt(2.0*Foam::constant::physicoChemical::k.value()/(rhoM_[cell]/rhoN_[cell])*translationalT_[cell]);

                instKnRho = MFP*maxMagGradRho[cell]/rhoM_[cell];
                instKnT = MFP*maxMagGradT[cell]/translationalT_[cell];
                instKnU = MFP*maxMagGradU[cell]/max(mag(UMean_[cell]),u0);

                instKnGLL = std::max(instKnRho,instKnT);
                instKnGLL = std::max(instKnGLL,instKnU); 

            }
            else
            {
                    instKnRho = 2.0*breakdownMax_;
                    instKnT = 2.0*breakdownMax_;
                    instKnU = 2.0*breakdownMax_;
                    instKnGLL = 2.0*breakdownMax_;
            }

            // time average macroscopic quantities
            KnRho_[cell] = theta_*instKnRho + (1.0-theta_)*KnRho_[cell];
            KnT_[cell] = theta_*instKnT + (1.0-theta_)*KnT_[cell];
            KnU_[cell] = theta_*instKnU + (1.0-theta_)*KnU_[cell];
            KnGLL_[cell] = theta_*instKnGLL + (1.0-theta_)*KnGLL_[cell];

        }

        KnRho_.correctBoundaryConditions();
        KnT_.correctBoundaryConditions();
        KnU_.correctBoundaryConditions();
        KnGLL_.correctBoundaryConditions();

        for (label pass=1; pass<=smoothingPasses_; pass++)
        {
            KnRho_ = fvc::average(fvc::interpolate(KnRho_));
            KnT_ = fvc::average(fvc::interpolate(KnT_));
            KnU_ = fvc::average(fvc::interpolate(KnU_));
            KnGLL_ = fvc::average(fvc::interpolate(KnGLL_)); 
            KnRho_.correctBoundaryConditions();
            KnT_.correctBoundaryConditions();
            KnU_.correctBoundaryConditions();
            KnGLL_.correctBoundaryConditions();
        }

        // determine cell collision model
        forAll(mesh_.cells(), cell)
        {
            if (KnGLL_[cell] > breakdownMax_)
            {
                cloud_.cellCollModel()[cell] = cloud_.binCollModel();
            }
            else
            {
                cloud_.cellCollModel()[cell] = cloud_.relCollModel();
            }
        }

        //Refine mesh decomposition
        label adjacentBinCollCells;
        label adjacentRelCollCells;
        label neighborRelCollCells;
        label neighborBinCollCells;

        for (label pass=1; pass<=refinementPasses_; pass++)
        {

            //Refine binary collision cells
            forAll(mesh_.cells(), cellI)
            {

                if (cloud_.cellCollModel()[cellI] == cloud_.binCollModel())
                {
                    adjacentBinCollCells = 0;
                    adjacentRelCollCells = 0;

                    forAll(mesh_.cellCells()[cellI], cellJ)
                    {
                        if (cloud_.cellCollModel()[mesh_.cellCells()[cellI][cellJ]] == cloud_.binCollModel())
                        {
                            adjacentBinCollCells++;
                        }
                        else
                        {
                            adjacentRelCollCells++;
                        }
                    }

                    if (adjacentBinCollCells == 0 || (adjacentBinCollCells == 1 && adjacentRelCollCells > 1))
                    {
                        cloud_.cellCollModel()[cellI] = cloud_.relCollModel();
                        continue;
                    }

                    fetchCellNeighborhood
                    (
                        cellI,
                        neighborLevels_,
                        neighborCells_
                    );

                    neighborBinCollCells = 0;
                    forAll(neighborCells_, cellJ)
                    {
                        if (cloud_.cellCollModel()[neighborCells_[cellJ]] == cloud_.binCollModel())
                        {
                            neighborBinCollCells++;
                        }
                    }  

                    if (neighborBinCollCells < maxNeighborFraction_*neighborCells_.size())
                    {
                        cloud_.cellCollModel()[cellI] = cloud_.relCollModel();
                        continue;
                    }

                }
            }

            //Refine relaxation collision cells
            forAll(mesh_.cells(), cellI)
            {

                if (cloud_.cellCollModel()[cellI] == cloud_.relCollModel())
                {
                
                    adjacentBinCollCells=0;
                    adjacentRelCollCells=0;

                    forAll(mesh_.cellCells()[cellI], cellJ)
                    {
                        if (cloud_.cellCollModel()[mesh_.cellCells()[cellI][cellJ]] == cloud_.binCollModel())
                        {
                            adjacentBinCollCells++;
                        }
                        else
                        {
                            adjacentRelCollCells++;
                        }
                    }

                    if (adjacentRelCollCells == 0 || (adjacentRelCollCells == 1 && adjacentBinCollCells > 1))
                    {
                        cloud_.cellCollModel()[cellI] = cloud_.binCollModel();
                        continue; 
                    }

                    fetchCellNeighborhood
                    (
                        cellI,
                        neighborLevels_,
                        neighborCells_
                    );

                    neighborRelCollCells = 0;
                    forAll(neighborCells_, cellJ)
                    {                    
                        if (cloud_.cellCollModel()[neighborCells_[cellJ]] == cloud_.relCollModel())
                        {
                            neighborRelCollCells++;
                        }
                    }  

                    if (neighborRelCollCells < maxNeighborFraction_*neighborCells_.size())
                    {
                        cloud_.cellCollModel()[cellI] = cloud_.binCollModel();
                        continue;
                    }

                }
            }
        }

        // reset
        timeSteps_ = 0;

        if (resetAtDecomposition_ && mesh_.time().value() < resetAtDecompositionUntilTime_+0.5*cloud_.mesh().time().deltaTValue())
        {

            timeAvCounter_ = 0.0;

            // reset cell information
            forAll(rhoN_, cell)
            {

                rhoNMean_[cell] = 0.0;
                rhoMMean_[cell] = 0.0;
                linearKEMean_[cell] = 0.0;
                rhoNMeanXnParticle_[cell] = 0.0;
                rhoMMeanXnParticle_[cell] = 0.0;
                linearKEMeanXnParticle_[cell] = 0.0;
                momentumMeanXnParticle_[cell] = vector::zero;

                forAll(typeIds_, i)
                {
                    nParcelsXnParticle_[i][cell] = 0.0;
                }
            }
            
        }

        update();
   
    }

}


void Foam::localKnudsen::update()
{

    // The main properties should be updated first
    updateProperties();

    propsDict_ = hybridDecompositionDict_.subDict(typeName + "Properties");

    propsDict_.readIfPresent("breakdownMax", breakdownMax_);

    propsDict_.readIfPresent("theta", theta_);

    propsDict_.readIfPresent("smoothingPasses", smoothingPasses_);  

}

// ************************************************************************* //

