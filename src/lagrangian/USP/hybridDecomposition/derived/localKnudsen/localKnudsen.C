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
    timeInterval_(dict.subDict("decompositionProperties").get<label>("timeInterval")),
    breakdownMax_(dict.subDict("decompositionProperties").get<scalar>("breakdownMax")),
    theta_(dict.subDict("decompositionProperties").getOrDefault<scalar>("theta",1.0)),
    smoothingPasses_(dict.subDict("decompositionProperties").getOrDefault<scalar>("smoothingPasses",0)),
    Tref_(dict.subDict("collisionProperties").get<scalar>("Tref")),
    timeSteps_(0),
    nAvTimeSteps_(0),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    linearKEMeanXnParticle_(mesh_.nCells(), 0.0),
    momentumMeanXnParticle_(mesh_.nCells(), Zero),
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
        dimensionedScalar(dimless, Zero)
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
        dimensionedScalar(dimless, Zero)
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
        dimensionedScalar(dimless, Zero)
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
        dimensionedScalar(dimless, Zero)
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
        dimensionedScalar(dimless/dimVolume, Zero)
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
        dimensionedScalar(dimMass/dimVolume, Zero)
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

    if (timeSteps_ == timeInterval_)
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
            rhoM_.correctBoundaryConditions();
            p_.correctBoundaryConditions();
            translationalT_.correctBoundaryConditions();
            UMean_.correctBoundaryConditions();
        }

        // assume single species - fix for mixtures
        const scalar& mass = cloud_.constProps(0).mass();
        const scalar& omega = cloud_.constProps(0).omega();
        const scalar& diameter = cloud_.constProps(0).d();
        const scalar& rotDoF = cloud_.constProps(0).rotationalDoF();
        scalar gamma = (5.0+rotDoF)/(3.0+rotDoF);

        scalarField magGradRho(mesh_.nCells());
        scalarField magGradT(mesh_.nCells());
        scalarField magGradU(mesh_.nCells());

        magGradRho = mag(fvc::grad(rhoM_));
        magGradT = mag(fvc::grad(translationalT_));
        magGradU = mag(fvc::grad(UMean_));

        // Calculate KnGLL based on macroscopic quantities
        scalar instKnRho;
        scalar instKnT;
        scalar instKnU;
        scalar instKnGLL;

        forAll(mesh_.cells(), cell)
        {

            if (rhoNMean_[cell] > VSMALL && translationalT_[cell] > VSMALL)
            {

                scalar meanFreePath = 1.0/(std::sqrt(2.0)*Foam::constant::mathematical::pi*sqr(diameter*std::pow(translationalT_[cell]/Tref_,0.5-omega))*(rhoM_[cell]/mass));
                scalar soundSpeed = std::sqrt(gamma*Foam::constant::physicoChemical::k.value()/mass*translationalT_[cell]);

                instKnRho = meanFreePath*magGradRho[cell]/rhoM_[cell];
                instKnT = meanFreePath*magGradT[cell]/translationalT_[cell];
                instKnU = meanFreePath*magGradU[cell]/max(mag(UMean_[cell]),soundSpeed);

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

