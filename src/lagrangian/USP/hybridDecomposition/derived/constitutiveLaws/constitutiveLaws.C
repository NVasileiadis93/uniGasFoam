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

#include "constitutiveLaws.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{

defineTypeNameAndDebug(constitutiveLaws, 0);

addToRunTimeSelectionTable(uspHybridDecomposition, constitutiveLaws, dictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveLaws::constitutiveLaws
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
    muu_(mesh_.nCells(), 0.0),
    muv_(mesh_.nCells(), 0.0),
    muw_(mesh_.nCells(), 0.0),
    mvv_(mesh_.nCells(), 0.0),
    mvw_(mesh_.nCells(), 0.0),
    mww_(mesh_.nCells(), 0.0),
    mcc_(mesh_.nCells(), 0.0),
    mccu_(mesh_.nCells(), 0.0),
    mccv_(mesh_.nCells(), 0.0),
    mccw_(mesh_.nCells(), 0.0),
    eu_(mesh_.nCells(), 0.0),
    ev_(mesh_.nCells(), 0.0),
    ew_(mesh_.nCells(), 0.0),
    e_(mesh_.nCells(), 0.0),
    nColls_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    linearKEMeanXnParticle_(mesh_.nCells(), 0.0),
    momentumMeanXnParticle_(mesh_.nCells(), Zero),
    nParcels_(),
    boundaryCells_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    momentumBF_(),
    CLB_
    (
        IOobject
        (
            "CLB",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    CLBQ_
    (
        IOobject
        (
            "CLBQ",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    CLBS_
    (
        IOobject
        (
            "CLBS",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    mu_
    (
        IOobject
        (
            "mu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure*dimTime, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass*dimLength/pow(dimTime,3)/dimTemperature, Zero),
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
    ),
    heatFluxVector_
    (
        IOobject
        (
            "heatFluxVector",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimMass*pow(dimTime,-3), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    pressureTensor_
    (
        IOobject
        (
            "pressureTensor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero)
    ),
    shearStressTensor_
    (
        IOobject
        (
            "shearStressTensor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    heatFluxVectorNS_
    (
        IOobject
        (
            "heatFluxVectorNS",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimMass*pow(dimTime,-3), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    shearStressTensorNS_
    (
        IOobject
        (
            "shearStressTensorNS",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero),
        zeroGradientFvPatchScalarField::typeName
    )
{

    typeIds_.setSize(cloud_.typeIdList().size());
    forAll(typeIds_, iD)
    {
        typeIds_[iD] = iD;
    }

    nParcels_.setSize(typeIds_.size());

    for (auto& n : nParcels_)
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

void Foam::constitutiveLaws::decompose()
{

 
    timeSteps_++;

    nAvTimeSteps_++;

    const scalar nParticle = cloud_.nParticle();

    // get cell measurements
    auto& cm = cloud_.cellPropMeasurements();

    forAll(typeIds_, iD)
    {

        rhoNMean_ += cm.rhoNMean()[iD];
        rhoMMean_  += cm.rhoMMean()[iD];
        linearKEMean_ += cm.linearKEMean()[iD];

        muu_ += cm.muu()[iD];
        muv_ += cm.muv()[iD];
        muw_ += cm.muw()[iD];
        mvv_ += cm.mvv()[iD];
        mvw_ += cm.mvw()[iD];
        mww_ += cm.mww()[iD];
        mcc_ += cm.mcc()[iD];
        mccu_ += cm.mccu()[iD];
        mccv_ += cm.mccv()[iD];
        mccw_ += cm.mccw()[iD];

        eu_ += cm.eu()[iD];
        ev_ += cm.ev()[iD];
        ew_ += cm.ew()[iD];
        e_ += cm.e()[iD];
        nColls_ += cm.nColls()[iD];

        rhoNMeanXnParticle_ += cm.rhoNMeanXnParticle()[iD];
        rhoMMeanXnParticle_ += cm.rhoMMeanXnParticle()[iD];
        linearKEMeanXnParticle_ += cm.linearKEMeanXnParticle()[iD];
        momentumMeanXnParticle_ += cm.momentumMeanXnParticle()[iD];

        nParcels_[iD] += cm.nParcels()[iD];

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

        scalar pSum;
        scalar tau;
        scalar Prandtl;

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

                // pressure tensor
                pressureTensor_[cell].xx() = rhoN_[cell]*
                (
                    muu_[cell]/(rhoNMean_[cell]) -
                    (
                        (rhoMMean_[cell]/(rhoNMean_[cell]))
                        *UMean_[cell].x()*UMean_[cell].x()
                    )
                );
                pressureTensor_[cell].xy() = rhoN_[cell]*
                (
                    muv_[cell]/(rhoNMean_[cell]) -
                    ((rhoMMean_[cell]/(rhoNMean_[cell])))
                    *UMean_[cell].x()*UMean_[cell].y()

                );
                pressureTensor_[cell].xz() = rhoN_[cell]*
                (
                    muw_[cell]/(rhoNMean_[cell]) -
                    ((rhoMMean_[cell]/(rhoNMean_[cell]))
                    *UMean_[cell].x()*UMean_[cell].z())
                );
                pressureTensor_[cell].yx() = pressureTensor_[cell].xy();
                pressureTensor_[cell].yy() = rhoN_[cell]*
                (
                    mvv_[cell]/(rhoNMean_[cell]) -
                    ((rhoMMean_[cell]/(rhoNMean_[cell])))
                    *UMean_[cell].y()*UMean_[cell].y()
                );
                pressureTensor_[cell].yz() = rhoN_[cell]*
                (
                    mvw_[cell]/(rhoNMean_[cell]) -
                    ((rhoMMean_[cell]/(rhoNMean_[cell]))
                    *UMean_[cell].y()*UMean_[cell].z())
                );
                pressureTensor_[cell].zx() = pressureTensor_[cell].xz();
                pressureTensor_[cell].zy() = pressureTensor_[cell].yz();
                pressureTensor_[cell].zz() = rhoN_[cell]*
                (
                    mww_[cell]/(rhoNMean_[cell]) -
                    ((rhoMMean_[cell]/(rhoNMean_[cell]))
                    *UMean_[cell].z()*UMean_[cell].z())
                );

                // Shear stress tensor
                scalar scalarPressure = (1.0/3.0)*
                                        (pressureTensor_[cell].xx() +
                                        pressureTensor_[cell].yy() +
                                        pressureTensor_[cell].zz());

                shearStressTensor_[cell] = -pressureTensor_[cell];
                shearStressTensor_[cell].xx() += scalarPressure;
                shearStressTensor_[cell].yy() += scalarPressure;
                shearStressTensor_[cell].zz() += scalarPressure;

                // Heat flux vector
                heatFluxVector_[cell].x() = rhoN_[cell]*
                (
                    0.5*(mccu_[cell]/(rhoNMean_[cell])) -
                    0.5*(mcc_[cell]/(rhoNMean_[cell]))*
                    UMean_[cell].x() + eu_[cell]/(rhoNMean_[cell]) -
                    (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].x()
                ) -
                    pressureTensor_[cell].xx()*UMean_[cell].x() -
                    pressureTensor_[cell].xy()*UMean_[cell].y() -
                    pressureTensor_[cell].xz()*UMean_[cell].z();

                heatFluxVector_[cell].y() = rhoN_[cell]*
                (
                    0.5*(mccv_[cell]/(rhoNMean_[cell])) -
                    0.5*(mcc_[cell]/(rhoNMean_[cell]))*
                    UMean_[cell].y() + ev_[cell]/(rhoNMean_[cell])-
                    (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].y()
                ) -
                    pressureTensor_[cell].yx()*UMean_[cell].x() -
                    pressureTensor_[cell].yy()*UMean_[cell].y() -
                    pressureTensor_[cell].yz()*UMean_[cell].z();

                heatFluxVector_[cell].z() = rhoN_[cell]*
                (
                    0.5*(mccw_[cell]/(rhoNMean_[cell])) -
                    0.5*(mcc_[cell]/(rhoNMean_[cell]))*
                    UMean_[cell].z() + ew_[cell]/(rhoNMean_[cell]) -
                    (e_[cell]/(rhoNMean_[cell]))*UMean_[cell].z()
                ) -
                    pressureTensor_[cell].zx()*UMean_[cell].x() -
                    pressureTensor_[cell].zy()*UMean_[cell].y() -
                    pressureTensor_[cell].zz()*UMean_[cell].z();

                //scale heat flux, pressure tensor and shear stress tensor for USP scheme
                pSum = 0.0;
                Prandtl = 0.0;
                forAll(typeIds_, iD)
                {
                    const scalar& rotDoF = cloud_.constProps(iD).rotationalDoF();
                    pSum += nParcels_[iD][cell];
                    Prandtl += nParcels_[iD][cell]*(5.0+rotDoF)/(7.5+rotDoF);
                }
                Prandtl /= pSum;
                tau = -0.5*log(1.0-pSum*nColls_[cell]/sqr(rhoNMean_[cell]));

                heatFluxVector_[cell] /= (1.0 - Prandtl*tau);
                pressureTensor_[cell] /= (1.0 - tau);
                shearStressTensor_[cell] /= (1.0 - tau);

            }
            else
            {
                rhoN_[cell] = 0.0;
                rhoM_[cell] = 0.0;
                p_[cell] = 0.0;
                translationalT_[cell] = 0.0;
                UMean_[cell] = vector::zero;
                heatFluxVector_[cell] = vector::zero;
                pressureTensor_[cell] = tensor::zero;
                shearStressTensor_[cell] = tensor::zero;
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

        // calculate heatflux and shear stress based on constitutive laws
        const scalar& mass = cloud_.constProps(0).mass();
        const scalar& omega = cloud_.constProps(0).omega();
        const scalar& diameter = cloud_.constProps(0).d();
        const scalar& rotDoF = cloud_.constProps(0).rotationalDoF();
        Prandtl = (5+rotDoF)/(7.5+rotDoF);
        const scalar CP = 2.5*Foam::constant::physicoChemical::k.value()/mass;
        const scalar muRef  = 7.5*std::sqrt(mass*Foam::constant::physicoChemical::k.value()*cloud_.collTref())
                       /(std::sqrt(Foam::constant::mathematical::pi)*(5.0-2.0*omega)*(7.0-2.0*omega)*sqr(diameter));

        forAll(mesh_.cells(), cell)
        {
            mu_[cell] = muRef*std::pow(translationalT_[cell]/cloud_.collTref(),omega);
            kappa_[cell] = mu_[cell]*CP/Prandtl;
        }

        volVectorField gradT = fvc::grad(translationalT_); 
        volTensorField gradU = fvc::grad(UMean_);
        gradT = fvc::grad(translationalT_);
        gradU = fvc::grad(UMean_);

        heatFluxVectorNS_ = -kappa_*gradT;
        shearStressTensorNS_ = mu_*(gradU + gradU.T()-(2e0/3e0)*I*tr(gradU));

        // correct boundary conditions
        heatFluxVector_.correctBoundaryConditions();
        heatFluxVectorNS_.correctBoundaryConditions();
        shearStressTensor_.correctBoundaryConditions();
        shearStressTensorNS_.correctBoundaryConditions();

        // smooth macroscopic quantities
        for (label pass=1; pass<=smoothingPasses_; pass++)
        {
            heatFluxVector_ = fvc::average(fvc::interpolate(heatFluxVector_));
            shearStressTensor_ = fvc::average(fvc::interpolate(shearStressTensor_));
            heatFluxVectorNS_ = fvc::average(fvc::interpolate(heatFluxVectorNS_));
            shearStressTensorNS_ = fvc::average(fvc::interpolate(shearStressTensorNS_));
            heatFluxVector_.correctBoundaryConditions();
            shearStressTensor_.correctBoundaryConditions();
            heatFluxVectorNS_.correctBoundaryConditions();
            shearStressTensorNS_.correctBoundaryConditions();                
        }

        scalar instCLB;
        vector instCLBQ;
        tensor instCLBS;
        forAll(mesh_.cells(), cell)
        {

            if (rhoNMean_[cell] > VSMALL && translationalT_[cell] > VSMALL)
            {
                scalar u0 = std::sqrt(2.0*Foam::constant::physicoChemical::k.value()/mass*translationalT_[cell]);

                instCLB = 0.0;
                forAll(instCLBQ, i) 
                {
                    instCLBQ[i] = fabs(heatFluxVector_[cell][i]-heatFluxVectorNS_[cell][i])/(p_[cell]*u0);
                    instCLB = max(instCLB, instCLBQ[i]);
                }

                forAll(instCLBS, i) 
                {
                    instCLBS[i] = fabs(shearStressTensor_[cell][i]-shearStressTensorNS_[cell][i])/(p_[cell]*u0);
                    instCLB = max(instCLB, instCLBS[i]);
                }
            }
            else
            {
                instCLB = GREAT;
                instCLBQ = vector(GREAT, GREAT, GREAT);
                instCLBS = tensor(GREAT, GREAT, GREAT, GREAT, GREAT, GREAT, GREAT, GREAT, GREAT);
            }

            // time average macroscopic quantities
            CLB_[cell] = theta_*instCLB + (1.0-theta_)*CLB_[cell];
            CLBQ_[cell] = theta_*instCLBQ + (1.0-theta_)*CLBQ_[cell];
            CLBS_[cell] = theta_*instCLBS + (1.0-theta_)*CLBS_[cell];

        }

        CLB_.correctBoundaryConditions();
        CLBQ_.correctBoundaryConditions();
        CLBS_.correctBoundaryConditions();

        // determine cell collision model
        forAll(mesh_.cells(), cell)
        {
            if (CLB_[cell] > breakdownMax_)
            {
                cloud_.cellCollModel()[cell] = cloud_.binCollModel();
            }
            else
            {
                cloud_.cellCollModel()[cell] = cloud_.relCollModel();
            }
        }

        //Remove isolated or single face connected cells
        for (label pass=1; pass<=10; pass++)
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
            muu_[cell] = 0.0;
            muv_[cell] = 0.0;
            muw_[cell] = 0.0;
            mvv_[cell] = 0.0;
            mvw_[cell] = 0.0;
            mww_[cell] = 0.0;
            mcc_[cell] = 0.0;
            mccu_[cell] = 0.0;
            mccv_[cell] = 0.0;
            mccw_[cell] = 0.0;
            eu_[cell] = 0.0;
            ev_[cell] = 0.0;
            ew_[cell] = 0.0;
            e_[cell] = 0.0;
            nColls_[cell] = 0.0;
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

