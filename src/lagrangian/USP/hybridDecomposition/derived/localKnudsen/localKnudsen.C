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
    const dictionary& dict,
    const polyMesh& mesh,
    uspCloud& cloud
)
:
    uspHybridDecomposition(dict, mesh, cloud),
    timeSteps_(0),
    nAvTimeSteps_(0),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    linearKEMeanXnParticle_(mesh_.nCells(), 0.0),
    momentumMeanXnParticle_(mesh_.nCells(), Zero),
    nParcelsXnParticle_(),
    boundaryCells_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    momentumBF_(),
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
    rhoM_
    (
        IOobject
        (
            "rhoM",
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
            "p",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure, Zero)
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
        dimensionedScalar(dimTemperature, Zero)
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
        dimensionedVector(dimVelocity, Zero)
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
        n.setSize(mesh_.nCells());
    }

    boundaryCells_.setSize(mesh_.boundaryMesh().size());

    forAll(boundaryCells_, p)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[p];

        boundaryCells_[p].setSize(patch.size());

        forAll(boundaryCells_[p], c)
        {
            boundaryCells_[p][c] = patch.faceCells()[c];
        }
    }

    // initialisation
    rhoNBF_.setSize(mesh_.boundaryMesh().size());
    rhoMBF_.setSize(mesh_.boundaryMesh().size());
    linearKEBF_.setSize(mesh_.boundaryMesh().size());
    momentumBF_.setSize(mesh_.boundaryMesh().size());

    forAll(rhoNBF_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];

        rhoNBF_[j].setSize(patch.size(), 0.0);
        rhoMBF_[j].setSize(patch.size(), 0.0);
        linearKEBF_[j].setSize(patch.size(), 0.0);
        momentumBF_[j].setSize(patch.size(), Zero);
    }

}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::localKnudsen::decompose()
{

    timeSteps_++;

    nAvTimeSteps_++;

    const scalar nParticle = cloud_.nParticle();

    // get cell measurements
    auto& cm = cloud_.cellPropMeasurements();

    forAll(cm.rhoNMean(), i)
    {

        rhoNMean_ += cm.rhoNMean()[i];
        rhoMMean_  += cm.rhoMMean()[i];
        linearKEMean_ += cm.linearKEMean()[i];
        rhoNMeanXnParticle_ += cm.rhoNMeanXnParticle()[i];
        rhoMMeanXnParticle_ += cm.rhoMMeanXnParticle()[i];
        linearKEMeanXnParticle_ += cm.linearKEMeanXnParticle()[i];
        momentumMeanXnParticle_ += cm.momentumMeanXnParticle()[i];
        nParcelsXnParticle_[i] += cm.nParcelsXnParticle()[i];

    }

    // Obtain boundary measurements
    auto& bm = cloud_.boundaryFluxMeasurements();

    forAll(bm.rhoNBF(), i)
    {

        forAll(bm.rhoNBF()[i], j)
        {
            forAll(bm.rhoNBF()[i][j], k)
            {
                rhoNBF_[j][k] += bm.rhoNBF()[i][j][k];
                rhoMBF_[j][k] += bm.rhoMBF()[i][j][k];
                linearKEBF_[j][k] += bm.linearKEBF()[i][j][k];
                momentumBF_[j][k] += bm.momentumBF()[i][j][k];
            }
        }
        
    }

    if (timeSteps_ == decompositionInterval_)
    {

        // computing internal fields
        forAll(rhoNMean_, cell)
        {
            if (rhoNMean_[cell] > VSMALL)
            {
                const scalar cellVolume = mesh_.cellVolumes()[cell];

                rhoN_[cell] = rhoNMeanXnParticle_[cell]/(nAvTimeSteps_*cellVolume);
                rhoM_[cell] = rhoMMeanXnParticle_[cell]/(nAvTimeSteps_*cellVolume);

                scalar rhoMMean = rhoMMeanXnParticle_[cell]/(cellVolume*nAvTimeSteps_);
                UMean_[cell] = momentumMeanXnParticle_[cell]/(rhoMMean*cellVolume*nAvTimeSteps_);

                scalar linearKEMean = 0.5*linearKEMeanXnParticle_[cell]/(cellVolume*nAvTimeSteps_);
                scalar rhoNMean = rhoNMeanXnParticle_[cell]/(cellVolume*nAvTimeSteps_);
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

        // Calcualte boundary vol fields
        forAll(rhoNBF_, j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];

            if (isA<wallPolyPatch>(patch))
            {
                forAll(rhoN_.boundaryFieldRef()[j], k)
                {

                    if (rhoNBF_[j][k] > VSMALL)
                    {

                        rhoN_.boundaryFieldRef()[j][k] = rhoNBF_[j][k]*nParticle/nAvTimeSteps_;
                        rhoM_.boundaryFieldRef()[j][k] = rhoMBF_[j][k]*nParticle/nAvTimeSteps_;
                        UMean_.boundaryFieldRef()[j][k] = momentumBF_[j][k]*nParticle/(rhoM_.boundaryFieldRef()[j][k]*nAvTimeSteps_);
                        scalar rhoMMean = rhoMBF_[j][k]*nParticle/nAvTimeSteps_;
                        scalar linearKEMean = linearKEBF_[j][k]*nParticle/nAvTimeSteps_;
                        scalar rhoNMean = rhoNBF_[j][k]*nParticle/nAvTimeSteps_;
                        translationalT_.boundaryFieldRef()[j][k] =
                            2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                            *(linearKEMean - 0.5*rhoMMean*(UMean_.boundaryFieldRef()[j][k] & UMean_.boundaryFieldRef()[j][k]));

                    }
                    else
                    {
                        rhoN_.boundaryFieldRef()[j][k] = 0.0;
                        rhoM_.boundaryFieldRef()[j][k] = 0.0;
                        translationalT_.boundaryFieldRef()[j][k] = 0.0;
                        UMean_.boundaryFieldRef()[j][k] = vector::zero;
                    }

                }

            }
        }

        forAll(boundaryCells_, j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];

            const labelList& bCs = boundaryCells_[j];

            forAll(bCs, k)
            {
                if
                (
                    isA<polyPatch>(patch)
                 && !isA<emptyPolyPatch>(patch)
                 && !isA<cyclicPolyPatch>(patch)
                )
                {

                    rhoN_.boundaryFieldRef()[j][k] = rhoN_[bCs[k]];
                    rhoM_.boundaryFieldRef()[j][k] = rhoM_[bCs[k]];

                    if (!isA<wallPolyPatch>(patch))
                    {
                        translationalT_.boundaryFieldRef()[j][k] =  translationalT_[bCs[k]];
                        p_.boundaryFieldRef()[j][k] = p_[bCs[k]];
                        UMean_.boundaryFieldRef()[j][k] = UMean_[bCs[k]];
                    }
                }
            }
        }

        // smooth macroscopic quantities
        for (label pass=1; pass<=smoothingPasses_; pass++)
        {

            rhoM_ = fvc::average(fvc::interpolate(rhoM_));
            p_ = fvc::average(fvc::interpolate(p_));
            translationalT_ = fvc::average(fvc::interpolate(translationalT_));
            UMean_ = fvc::average(fvc::interpolate(UMean_)); 

            forAll(boundaryCells_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];

                const labelList& bCs = boundaryCells_[j];

                forAll(bCs, k)
                {
                    if
                    (
                        isA<polyPatch>(patch)
                     && !isA<emptyPolyPatch>(patch)
                     && !isA<cyclicPolyPatch>(patch)
                    )
                    {
                        rhoM_.boundaryFieldRef()[j][k] = boundCoeff_*rhoM_.boundaryFieldRef()[j][k]+(1.0-boundCoeff_)*rhoM_[bCs[k]];
                        translationalT_.boundaryFieldRef()[j][k] =  boundCoeff_*translationalT_.boundaryFieldRef()[j][k]+(1.0-boundCoeff_)*translationalT_[bCs[k]];
                        p_.boundaryFieldRef()[j][k] = boundCoeff_*p_.boundaryFieldRef()[j][k]+(1.0-boundCoeff_)*p_[bCs[k]];
                        UMean_.boundaryFieldRef()[j][k] = boundCoeff_*UMean_.boundaryFieldRef()[j][k]+(1.0-boundCoeff_)*UMean_[bCs[k]];
                    }
                }
            }

        }

        // assume single species - fix for mixtures
        const scalar& mass = cloud_.constProps(0).mass();
        const scalar& omega = cloud_.constProps(0).omega();
        const scalar& diameter = cloud_.constProps(0).d();
        const scalar& rotDoF = cloud_.constProps(0).rotationalDoF();
        scalar gamma = (5.0+rotDoF)/(3.0+rotDoF);

        scalarField maxMagGradRho(mesh_.nCells());
        scalarField maxMagGradT(mesh_.nCells());
        scalarField maxMagGradU(mesh_.nCells());

        maxMagGradRho = mag(fvc::grad(rhoM_));
        maxMagGradT = mag(fvc::grad(translationalT_));
        maxMagGradU = mag(fvc::grad(UMean_));

        //forAll(mesh_.cells(), cell)
        //{
        //
        //    maxMagGradRho[cell] = 0.0;
        //    maxMagGradT[cell] = 0.0;
        //    maxMagGradU[cell] = 0.0;
        //    forAll(mesh_.cellCells()[cell], cellI)
        //    {
        //
        //        label adjCell = mesh_.cellCells()[cell][cellI];
        //        scalar cellDistance = mag(mesh_.C()[adjCell]-mesh_.C()[cell]);
        //
        //        maxMagGradRho[cell] = max(maxMagGradRho[cell], fabs(rhoM_[adjCell]-rhoM_[cell])/cellDistance);
        //        maxMagGradT[cell] = max(maxMagGradT[cell], fabs(translationalT_[adjCell]-translationalT_[cell])/cellDistance);
        //        maxMagGradU[cell] = max(maxMagGradU[cell], fabs(mag(UMean_[adjCell])-mag(UMean_[cell]))/cellDistance);
        //
        //    }
        //}        

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
                    label qspec = 0;

                    for (qspec=0; qspec<typeIds_.size(); ++qspec)
                    {
                        scalar dPQ = 0.5*(cloud_.constProps(i).d() + cloud_.constProps(qspec).d());

                        scalar omegaPQ = 0.5*(cloud_.constProps(i).omega() + cloud_.constProps(qspec).omega());

                        scalar massRatio = cloud_.constProps(i).mass()/cloud_.constProps(qspec).mass();

                        if (nParcelsXnParticle_[qspec][cell] > VSMALL)
                        {

                            scalar nDensQ = nParcelsXnParticle_[qspec][cell]/cellVolume;

                            scalar reducedMass =
                                cloud_.constProps(i).mass()*cloud_.constProps(qspec).mass()
                                /(cloud_.constProps(i).mass() + cloud_.constProps(qspec).mass());

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

                scalar u0 = std::sqrt(gamma*Foam::constant::physicoChemical::k.value()/(rhoM_[cell]/rhoN_[cell])*translationalT_[cell]);

                instKnRho = MFP*maxMagGradRho[cell]/rhoM_[cell];
                instKnT = MFP*maxMagGradT[cell]/translationalT_[cell];
                instKnU = MFP*maxMagGradU[cell]/max(mag(UMean_[cell]),u0);

                instKnGLL = std::max(instKnRho,instKnT);
                instKnGLL = std::max(instKnGLL,instKnU); 

            }
            else
            {
                    instKnRho = GREAT;
                    instKnT = GREAT;
                    instKnU = GREAT;
                    instKnGLL = GREAT;
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

        //Remove isolated or single face connected cells
        for (label pass=1; pass<=refinementPasses_; pass++)
        {
            forAll(mesh_.cells(), cell)
            {

                if (cloud_.cellCollModel()[cell] == cloud_.binCollModel())
                {

                    label adjacentBinCollCells=0;
                    
                    forAll(mesh_.cellCells()[cell], adjCell)
                    {
                        if (cloud_.cellCollModel()[mesh_.cellCells()[cell][adjCell]] == cloud_.binCollModel())
                        {
                            adjacentBinCollCells++;
                        }
                    }
                    
                    if (adjacentBinCollCells <= 1)
                    {
                        cloud_.cellCollModel()[cell] = cloud_.relCollModel();   
                    }

                }

                if (cloud_.cellCollModel()[cell] == cloud_.relCollModel())
                {

                    label adjacentRelCollCells=0;
                    
                    forAll(mesh_.cellCells()[cell], adjCell)
                    {
                        if (cloud_.cellCollModel()[mesh_.cellCells()[cell][adjCell]] == cloud_.relCollModel())
                        {
                            adjacentRelCollCells++;
                        }
                    }
                    
                    if (adjacentRelCollCells <= 1)
                    {
                        cloud_.cellCollModel()[cell] = cloud_.binCollModel();   
                    }

                }

            }
        }

        // reset
        timeSteps_ = 0;

        nAvTimeSteps_ = 0;

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

        // reset boundary information
        forAll(rhoNBF_, j)
        {
            rhoNBF_[j] = 0.0;
            rhoMBF_[j] = 0.0;
            linearKEBF_[j] = 0.0;
            momentumBF_[j] = vector::zero;
        }

    }

}

// ************************************************************************* //

