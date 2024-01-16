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

#include "uspHybridDecomposer.H"
#include "uspCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
uspHybridDecomposer::uspHybridDecomposer
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
    firstDecomp(true),
    timeSteps_(0),
    nAvTimeSteps_(0),
    breakdownCriterion_(dict.subDict("hybridProperties").get<word>("breakdownCriterion")),
    decomposeInterval_(dict.subDict("hybridProperties").get<label>("decomposeInterval")),
    breakdownMax_(dict.subDict("hybridProperties").get<scalar>("breakdownMax")),
    smoothingPasses_(dict.subDict("hybridProperties").get<scalar>("smoothingPasses")),
    refinementPasses_(dict.subDict("hybridProperties").get<scalar>("refinementPasses")),
    Tref_(dict.subDict("collisionCoeffs").get<scalar>("Tref")),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    linearKEMeanXnParticle_(mesh_.nCells(), 0.0),
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
    momentumMean_(mesh_.nCells(), Zero),
    momentumMeanXnParticle_(mesh_.nCells(), Zero),
    nParcels_(),
    nParcelsXnParticle_(),
    mccSpecies_(),
    boundaryCells_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    momentumBF_(),
    speciesRhoNBF_(),
    B_
    (
        IOobject
        (
            "B_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    KnGLL_
    (
        IOobject
        (
            "KnGLL_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    KnRho_
    (
        IOobject
        (
            "KnRho_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    KnT_
    (
        IOobject
        (
            "KnT_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    KnU_
    (
        IOobject
        (
            "KnU_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    CLB_
    (
        IOobject
        (
            "CLB",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
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
            IOobject::NO_READ,
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
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    ASq_
    (
        IOobject
        (
            "ASq_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    breakdownParameter_
    (
        IOobject
        (
            "breakdownParameter_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    mu_
    (
        IOobject
        (
            "mu",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure*dimTime, Zero)
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimMass*dimLength/pow(dimTime,3)/dimTemperature, Zero)
    ),
    prevRhoN_
    (
        IOobject
        (
            "prevRhoN_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, Zero)
    ),
    prevRhoM_
    (
        IOobject
        (
            "prevRhoM_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimVolume, Zero)
    ),
    prevP_
    (
        IOobject
        (
            "prevP_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure, Zero)
    ),
    prevTranslationalT_
    (
        IOobject
        (
            "prevTranslationalT_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, Zero)
    ),
    prevUMean_
    (
        IOobject
        (
            "prevUMean_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimVelocity, Zero)
    ),
    prevHeatFluxVector_
    (
        IOobject
        (
            "prevHeatFluxVector_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimMass*pow(dimTime,-3), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    prevPressureTensor_
    (
        IOobject
        (
            "prevPressureTensor_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero)
    ),
    prevShearStressTensor_
    (
        IOobject
        (
            "prevShearStressTensor_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero),
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
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, Zero)
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimVolume, Zero)
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
        dimensionedScalar(dimPressure, Zero)
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
        dimensionedScalar(dimTemperature, Zero)
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
        dimensionedVector(dimVelocity, Zero)
    ),
    heatFluxVector_
    (
        IOobject
        (
            "heatFluxVector_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimMass*pow(dimTime,-3), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    pressureTensor_
    (
        IOobject
        (
            "pressureTensor_",
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
            "shearStressTensor_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    heatFluxVectorNS_
    (
        IOobject
        (
            "heatFluxVectorNS_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimMass*pow(dimTime,-3), Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    shearStressTensorNS_
    (
        IOobject
        (
            "shearStressTensorNS_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero),
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

    boundaryCells_.setSize(mesh.boundaryMesh().size());

    forAll(boundaryCells_, p)
    {
        const polyPatch& patch = mesh.boundaryMesh()[p];

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

void uspHybridDecomposer::update()
{

    scalar theta_ = 0.2;

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
        momentumMean_ += cm.momentumMean()[i];
        rhoNMeanXnParticle_ += cm.rhoNMeanXnParticle()[i];
        rhoMMeanXnParticle_ += cm.rhoMMeanXnParticle()[i];
        momentumMeanXnParticle_ += cm.momentumMeanXnParticle()[i];
        linearKEMeanXnParticle_ += cm.linearKEMeanXnParticle()[i];

        muu_ += cm.muu()[i];
        muv_ += cm.muv()[i];
        muw_ += cm.muw()[i];
        mvv_ += cm.mvv()[i];
        mvw_ += cm.mvw()[i];
        mww_ += cm.mww()[i];
        mcc_ += cm.mcc()[i];
        mccu_ += cm.mccu()[i];
        mccv_ += cm.mccv()[i];
        mccw_ += cm.mccw()[i];

        eu_ += cm.eu()[i];
        ev_ += cm.ev()[i];
        ew_ += cm.ew()[i];
        e_ += cm.e()[i];

        nParcels_[i] += cm.nParcels()[i];
        nParcelsXnParticle_[i] += cm.nParcelsXnParticle()[i];

    }

    // Obtain boundary measurements
    auto& bm = cloud_.boundaryFluxMeasurements();

    forAll(bm.rhoNBF(), i)
    {
        const label iD = typeIds_.find(i);

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

    /*const label maxParticles = 20;
    List<DynamicList<uspParcel*>>& cellOccupancy = cloud_.cellOccupancy();

    forAll(mesh_.cells(), cellI)
    {
        //Anderson-Darling Test
        const DynamicList<uspParcel*>& cellParcels(cellOccupancy[cellI]);
        const label n = min(cellParcels.size(),maxParticles);

        scalar maxAD = 0.0;
        if (n > 1)
        {
            scalar instAD = 0.0;
            
            for (label j=0; j<3; j++) //// CHANGE j<3
            {

                SortableList<scalar> vList(n);

                forAll(vList, i)
                {
                    uspParcel& p = *cellParcels[i];
                    vList[i] = p.U()[j];
                }
                vList.sort();

                scalar vListAvg = 0.0;
                forAll(vList, i)
                {
                    vListAvg += vList[i];
                }
                vListAvg /= scalar(n);

                scalar vListStDev = 0.0;
                forAll(vList, i)
                {
                    vListStDev += sqr(vList[i]-vListAvg);
                }
                vListStDev = sqrt(vListStDev/scalar(n-1));


                if (vListStDev > 0)
                {
                    forAll(vList, i)
                    {
                        vList[i] = (vList[i]-vListAvg)/vListStDev;
                    }

                    instAD = 0.0;
                    scalar cumDistFunct;
                    forAll(vList, i)
                    {
                        cumDistFunct = 0.5*(1.0+erf(vList[i]/sqrt(2.0)));
                        //Info << vList[i] << " " << cumDistFunct << endl;
                        instAD += (2.0*i+1.0)*log(cumDistFunct)+(2.0*(n-i)-1.0)*log(1.0-cumDistFunct);
                        //instAD += sqr(0.5*(2.0*i+1.0)/scalar(n)-cumDistFunct);
                    }     
                    instAD = -n - instAD/scalar(n);
                    //instAD += 1.0/(12.0*scalar(n)); 
                    //Info << instAD << endl;
                    //std::cin.get();
                }
                else
                {
                    instAD += 10.0;
                }

                if (maxAD < instAD)
                {
                    maxAD = instAD;
                }

            }

        }
        else
        {
            maxAD += 10.0;
        }

        ASq_[cellI] += maxAD*(1.0+0.75/scalar(n)+2.25/sqr(scalar(n)));
        
    }*/

    if (timeSteps_ == decomposeInterval_)
    {

        // computing internal fields
        forAll(rhoNMean_, cellI)
        {
            if (rhoNMean_[cellI] > VSMALL)
            {
                const scalar cellVolume = mesh_.cellVolumes()[cellI];

                rhoN_[cellI] = rhoNMeanXnParticle_[cellI]/(nAvTimeSteps_*cellVolume);
                rhoM_[cellI] = rhoMMeanXnParticle_[cellI]/(nAvTimeSteps_*cellVolume);

                scalar rhoMMean = rhoMMeanXnParticle_[cellI]/(cellVolume*nAvTimeSteps_);
                UMean_[cellI] = momentumMeanXnParticle_[cellI]/(rhoMMean*cellVolume*nAvTimeSteps_);

                scalar linearKEMean = 0.5*linearKEMeanXnParticle_[cellI]/(cellVolume*nAvTimeSteps_);
                scalar rhoNMean = rhoNMeanXnParticle_[cellI]/(cellVolume*nAvTimeSteps_);
                translationalT_[cellI] =
                    2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                   *(
                        linearKEMean
                      - 0.5*rhoMMean*(UMean_[cellI] & UMean_[cellI])
                    );

                p_[cellI] = rhoN_[cellI]*physicoChemical::k.value()*translationalT_[cellI];

                // pressure tensor
                pressureTensor_[cellI].xx() = rhoN_[cellI]*
                (
                    muu_[cellI]/(rhoNMean_[cellI]) -
                    (
                        (rhoMMean_[cellI]/(rhoNMean_[cellI]))
                        *UMean_[cellI].x()*UMean_[cellI].x()
                    )
                );
                pressureTensor_[cellI].xy() = rhoN_[cellI]*
                (
                    muv_[cellI]/(rhoNMean_[cellI]) -
                    ((rhoMMean_[cellI]/(rhoNMean_[cellI])))
                    *UMean_[cellI].x()*UMean_[cellI].y()

                );
                pressureTensor_[cellI].xz() = rhoN_[cellI]*
                (
                    muw_[cellI]/(rhoNMean_[cellI]) -
                    ((rhoMMean_[cellI]/(rhoNMean_[cellI]))
                    *UMean_[cellI].x()*UMean_[cellI].z())
                );
                pressureTensor_[cellI].yx() = pressureTensor_[cellI].xy();
                pressureTensor_[cellI].yy() = rhoN_[cellI]*
                (
                    mvv_[cellI]/(rhoNMean_[cellI]) -
                    ((rhoMMean_[cellI]/(rhoNMean_[cellI])))
                    *UMean_[cellI].y()*UMean_[cellI].y()
                );
                pressureTensor_[cellI].yz() = rhoN_[cellI]*
                (
                    mvw_[cellI]/(rhoNMean_[cellI]) -
                    ((rhoMMean_[cellI]/(rhoNMean_[cellI]))
                    *UMean_[cellI].y()*UMean_[cellI].z())
                );
                pressureTensor_[cellI].zx() = pressureTensor_[cellI].xz();
                pressureTensor_[cellI].zy() = pressureTensor_[cellI].yz();
                pressureTensor_[cellI].zz() = rhoN_[cellI]*
                (
                    mww_[cellI]/(rhoNMean_[cellI]) -
                    ((rhoMMean_[cellI]/(rhoNMean_[cellI]))
                    *UMean_[cellI].z()*UMean_[cellI].z())
                );

                // Shear stress tensor
                scalar scalarPressure = (1.0/3.0)*
                                        (pressureTensor_[cellI].xx() +
                                        pressureTensor_[cellI].yy() +
                                        pressureTensor_[cellI].zz());

                shearStressTensor_[cellI] = -pressureTensor_[cellI];
                shearStressTensor_[cellI].xx() += scalarPressure;
                shearStressTensor_[cellI].yy() += scalarPressure;
                shearStressTensor_[cellI].zz() += scalarPressure;

                // Heat flux vector
                heatFluxVector_[cellI].x() = rhoN_[cellI]*
                (
                    0.5*(mccu_[cellI]/(rhoNMean_[cellI])) -
                    0.5*(mcc_[cellI]/(rhoNMean_[cellI]))*
                    UMean_[cellI].x() + eu_[cellI]/(rhoNMean_[cellI]) -
                    (e_[cellI]/(rhoNMean_[cellI]))*UMean_[cellI].x()
                ) -
                    pressureTensor_[cellI].xx()*UMean_[cellI].x() -
                    pressureTensor_[cellI].xy()*UMean_[cellI].y() -
                    pressureTensor_[cellI].xz()*UMean_[cellI].z();

                heatFluxVector_[cellI].y() = rhoN_[cellI]*
                (
                    0.5*(mccv_[cellI]/(rhoNMean_[cellI])) -
                    0.5*(mcc_[cellI]/(rhoNMean_[cellI]))*
                    UMean_[cellI].y() + ev_[cellI]/(rhoNMean_[cellI])-
                    (e_[cellI]/(rhoNMean_[cellI]))*UMean_[cellI].y()
                ) -
                    pressureTensor_[cellI].yx()*UMean_[cellI].x() -
                    pressureTensor_[cellI].yy()*UMean_[cellI].y() -
                    pressureTensor_[cellI].yz()*UMean_[cellI].z();

                heatFluxVector_[cellI].z() = rhoN_[cellI]*
                (
                    0.5*(mccw_[cellI]/(rhoNMean_[cellI])) -
                    0.5*(mcc_[cellI]/(rhoNMean_[cellI]))*
                    UMean_[cellI].z() + ew_[cellI]/(rhoNMean_[cellI]) -
                    (e_[cellI]/(rhoNMean_[cellI]))*UMean_[cellI].z()
                ) -
                    pressureTensor_[cellI].zx()*UMean_[cellI].x() -
                    pressureTensor_[cellI].zy()*UMean_[cellI].y() -
                    pressureTensor_[cellI].zz()*UMean_[cellI].z();

            }
            else
            {
                rhoN_[cellI] = 0.0;
                rhoM_[cellI] = 0.0;
                p_[cellI] = 0.0;
                translationalT_[cellI] = 0.0;
                UMean_[cellI] = vector::zero;
                heatFluxVector_[cellI] = vector::zero;
                pressureTensor_[cellI] = tensor::zero;
                shearStressTensor_[cellI] = tensor::zero;
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
                    rhoN_.boundaryFieldRef()[j][k] = rhoNBF_[j][k]*nParticle/nAvTimeSteps_;
                    rhoM_.boundaryFieldRef()[j][k] = rhoMBF_[j][k]*nParticle/nAvTimeSteps_;

                    if (rhoM_.boundaryFieldRef()[j][k] > VSMALL)
                    {
                        UMean_.boundaryFieldRef()[j][k] =
                            momentumBF_[j][k]*nParticle
                        /(rhoM_.boundaryFieldRef()[j][k]*nAvTimeSteps_);
                    }
                    else
                    {
                        UMean_.boundaryFieldRef()[j][k] = vector::zero;
                    }

                    scalar rhoMMean = rhoMBF_[j][k]*nParticle/nAvTimeSteps_;
                    scalar linearKEMean = linearKEBF_[j][k]*nParticle/nAvTimeSteps_;
                    scalar rhoNMean = rhoNBF_[j][k]*nParticle/nAvTimeSteps_;

                    if (rhoNMean > VSMALL)
                    {
                        translationalT_.boundaryFieldRef()[j][k] =
                            2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                        *(
                                linearKEMean
                            - 0.5*rhoMMean
                            *(
                                    UMean_.boundaryFieldRef()[j][k]
                                & UMean_.boundaryFieldRef()[j][k]
                                )
                            );
                    }
                    else
                    {
                        translationalT_.boundaryFieldRef()[j][k] = 0.0;
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


        if (!firstDecomp)
        {
            rhoN_ = theta_*rhoN_ + (1.0-theta_)*prevRhoN_;
            rhoM_ = theta_*rhoM_ + (1.0-theta_)*prevRhoM_;
            p_ = theta_*p_ + (1.0-theta_)*prevP_;
            translationalT_ = theta_*translationalT_ + (1.0-theta_)*prevTranslationalT_;
            UMean_ = theta_*UMean_ + (1.0-theta_)*prevUMean_ ;
            heatFluxVector_ = theta_*heatFluxVector_ + (1.0-theta_)*prevHeatFluxVector_;
            pressureTensor_ = theta_*pressureTensor_ + (1.0-theta_)*prevPressureTensor_;
            shearStressTensor_ = theta_*shearStressTensor_ + (1.0-theta_)*prevShearStressTensor_;
        }
        else
        {
            firstDecomp = false;
        }

        prevRhoN_ = rhoN_;
        prevRhoM_ = rhoM_;
        prevP_ = p_;
        prevTranslationalT_ = translationalT_;
        prevUMean_ = UMean_;
        prevHeatFluxVector_ = heatFluxVector_;
        prevPressureTensor_ = pressureTensor_;
        prevShearStressTensor_ = shearStressTensor_;

        const scalar& mass = cloud_.constProps(0).mass();
        const scalar& omega = cloud_.constProps(0).omega();
        const scalar& diameter = cloud_.constProps(0).d();
        const scalar& rotDoF = cloud_.constProps(0).rotationalDoF();
        const scalar Prandtl = (5+rotDoF)/(7.5+rotDoF);
        const scalar CP = 2.5*Foam::constant::physicoChemical::k.value()/mass;
        const scalar muRef  = 7.5*std::sqrt(mass*Foam::constant::physicoChemical::k.value()*Tref_)
                       /(std::sqrt(Foam::constant::mathematical::pi)*(5.0-2.0*omega)*(7.0-2.0*omega)*sqr(diameter));

        forAll(mesh_.cells(), cellI)
        {
            mu_[cellI] = muRef*std::pow(translationalT_[cellI]/Tref_,omega);
            kappa_[cellI] = mu_[cellI]*CP/Prandtl;
        }

        volVectorField gradT = fvc::grad(translationalT_); 
        volTensorField gradU = fvc::grad(UMean_);
        gradT = fvc::grad(translationalT_);
        gradU = fvc::grad(UMean_);

        heatFluxVectorNS_ = -kappa_*gradT;
        shearStressTensorNS_ = mu_*(gradU + gradU.T()-(2e0/3e0)*I*tr(gradU));

        if (breakdownCriterion_ == "ChapmanEnskog")
        {

            // smooth macroscopic quantities
            for (label pass=1; pass<=smoothingPasses_; pass++)
            {
                p_ = fvc::average(fvc::interpolate(p_));
                translationalT_ = fvc::average(fvc::interpolate(translationalT_));
                heatFluxVector_ = fvc::average(fvc::interpolate(heatFluxVector_));
                shearStressTensor_ = fvc::average(fvc::interpolate(shearStressTensor_));
                p_.correctBoundaryConditions();
                translationalT_.correctBoundaryConditions();
                heatFluxVector_.correctBoundaryConditions();
                shearStressTensor_.correctBoundaryConditions();
            }

            // Calculate breakdown parameter
            forAll(mesh_.cells(), cellI)
            {

                scalar massMean = 0.0;
                forAll(typeIds_,iD)
                {
                    massMean += cloud_.constProps(iD).mass()*nParcels_[iD][cellI];
                }
                massMean /= rhoNMean_[cellI];

                scalar mostProbSpeed(
                    cloud_.maxwellianMostProbableSpeed
                    (
                        translationalT_[cellI],
                        massMean
                    )
                );

                B_[cellI] = max(mag(heatFluxVector_[cellI])/mostProbSpeed,mag(shearStressTensor_[cellI]))/p_[cellI];

            }

            breakdownParameter_ = B_;

        }
        else if (breakdownCriterion_ == "localKnudsen")
        {

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

            // fix for mixtures
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
 
            // Calculate KnMax based on macroscopic quantities
            forAll(mesh_.cells(), cellI)
            {

                scalar meanFreePath = 1.0/(std::sqrt(2.0)*Foam::constant::mathematical::pi*sqr(diameter*std::pow(translationalT_[cellI]/Tref_,0.5-omega))*(rhoM_[cellI]/mass));
                scalar soundSpeed = std::sqrt(gamma*Foam::constant::physicoChemical::k.value()/mass*translationalT_[cellI]);

                KnRho_[cellI] = meanFreePath*magGradRho[cellI]/rhoM_[cellI];
                KnT_[cellI] = meanFreePath*magGradT[cellI]/translationalT_[cellI];
                KnU_[cellI] = meanFreePath*magGradU[cellI]/max(mag(UMean_[cellI]),soundSpeed);

                KnGLL_[cellI]=std::max(KnRho_[cellI],KnT_[cellI]);
                KnGLL_[cellI]=std::max(KnGLL_[cellI],KnU_[cellI]); 

            }

            breakdownParameter_ = KnGLL_;

        } 
        else if (breakdownCriterion_ == "constitutiveLaw")
        {

            // normalise heatflux and stress
            scalar u0;
            forAll(mesh_.cells(), cellI)
            {

                u0 = std::sqrt(2.0*Foam::constant::physicoChemical::k.value()/mass*translationalT_[cellI]);

                heatFluxVector_[cellI] = heatFluxVector_[cellI]/(p_[cellI]*u0);
                heatFluxVectorNS_[cellI] = heatFluxVectorNS_[cellI]/(p_[cellI]*u0);
                shearStressTensor_[cellI] = shearStressTensor_[cellI]/p_[cellI];
                shearStressTensorNS_[cellI] = shearStressTensorNS_[cellI]/p_[cellI];

            }
            heatFluxVector_.correctBoundaryConditions();
            heatFluxVectorNS_.correctBoundaryConditions();
            shearStressTensor_.correctBoundaryConditions();
            shearStressTensorNS_.correctBoundaryConditions();

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

            const scalar bMax = 0.05;
            forAll(mesh_.cells(), cellI)
            {

                CLB_[cellI] = 0.0;
                forAll(CLBQ_[cellI],i) 
                {
                    CLBQ_[cellI][i] = sqr(heatFluxVector_[cellI][i]-heatFluxVectorNS_[cellI][i]);
                    CLB_[cellI] += CLBQ_[cellI][i]; 
                }

                forAll(CLBS_[cellI],i) 
                {
                    CLBS_[cellI][i] = sqr(shearStressTensor_[cellI][i]-shearStressTensorNS_[cellI][i]);
                    CLB_[cellI] += CLBS_[cellI][i];
                }

                CLB_[cellI] = std::sqrt(CLB_[cellI]/12.0);

            }

            breakdownParameter_ = CLB_;

        }
        else if (breakdownCriterion_ == "AndersonDarling")
        {

            breakdownParameter_ = ASq_/nAvTimeSteps_;

        }
        else
        {
            FatalErrorInFunction
                << "Continuum breakdown criterion "
                << breakdownCriterion_
                << " not found."
                << exit(FatalError);
        }

        // determine cell collision model
        forAll(mesh_.cells(), cellI)
        {
            if (breakdownParameter_[cellI] > breakdownMax_)
            {
                cloud_.cellCollModel()[cellI] = cloud_.binCollModel();
            }
            else
            {
                cloud_.cellCollModel()[cellI] = cloud_.relCollModel();
            }
        }

        //Remove isolated or single face connected cells
        for (label pass=1; pass<=100; pass++)
        {
            forAll(mesh_.cells(), cellI)
            {

                if (cloud_.cellCollModel()[cellI] == cloud_.binCollModel())
                {

                    label adjacentBinCollCells=0;
                    
                    forAll(mesh_.cellCells()[cellI], cellJ)
                    {
                        if (cloud_.cellCollModel()[mesh_.cellCells()[cellI][cellJ]] == cloud_.binCollModel())
                        {
                            adjacentBinCollCells++;
                        }
                    }
                    
                    if (adjacentBinCollCells <= 1)
                    {
                        cloud_.cellCollModel()[cellI] = cloud_.relCollModel();   
                    }

                }

                if (cloud_.cellCollModel()[cellI] == cloud_.relCollModel())
                {

                    label adjacentRelCollCells=0;
                    
                    forAll(mesh_.cellCells()[cellI], cellJ)
                    {
                        if (cloud_.cellCollModel()[mesh_.cellCells()[cellI][cellJ]] == cloud_.relCollModel())
                        {
                            adjacentRelCollCells++;
                        }
                    }
                    
                    if (adjacentRelCollCells <= 1)
                    {
                        cloud_.cellCollModel()[cellI] = cloud_.binCollModel();   
                    }

                }

            }
        }

        // reset
        timeSteps_ = 0;

        nAvTimeSteps_ = 0;

        forAll(rhoN_, cellI)
        {

            rhoNMean_[cellI] = scalar(0.0);
            rhoMMean_[cellI] = scalar(0.0);
            linearKEMean_[cellI] = scalar(0.0);
            momentumMean_[cellI] = vector::zero;
            muu_[cellI] = scalar(0.0);
            muv_[cellI] = scalar(0.0);
            muw_[cellI] = scalar(0.0);
            mvv_[cellI] = scalar(0.0);
            mvw_[cellI] = scalar(0.0);
            mww_[cellI] = scalar(0.0);
            mcc_[cellI] = scalar(0.0);
            mccu_[cellI] = scalar(0.0);
            mccv_[cellI] = scalar(0.0);
            mccw_[cellI] = scalar(0.0);
            eu_[cellI] = scalar(0.0);
            ev_[cellI] = scalar(0.0);
            ew_[cellI] = scalar(0.0);
            e_[cellI] = scalar(0.0);
            rhoNMeanXnParticle_[cellI] = scalar(0.0);
            rhoMMeanXnParticle_[cellI] = scalar(0.0);
            momentumMeanXnParticle_[cellI] = vector::zero;
            linearKEMeanXnParticle_[cellI] = scalar(0.0);

            ASq_[cellI] = scalar(0.0);

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

}  // End namespace Foam

// ************************************************************************* //

