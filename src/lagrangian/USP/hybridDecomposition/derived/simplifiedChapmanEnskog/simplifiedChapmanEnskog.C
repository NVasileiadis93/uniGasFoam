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

#include "simplifiedChapmanEnskog.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{

defineTypeNameAndDebug(simplifiedChapmanEnskog, 0);

addToRunTimeSelectionTable(uspHybridDecomposition, simplifiedChapmanEnskog, dictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simplifiedChapmanEnskog::simplifiedChapmanEnskog
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
    B_
    (
        IOobject
        (
            "B",
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
        dimensionedScalar(dimPressure, Zero),
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
    heatFluxVector_
    (
        IOobject
        (
            "heatFluxVectorCE",
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
            "pressureTensor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    shearStressTensor_
    (
        IOobject
        (
            "shearStressTensorCE",
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
    forAll(typeIds_, iD)
    {
        typeIds_[iD] = iD;
    }

    nParcels_.setSize(typeIds_.size());

    for (auto& n : nParcels_)
    {
        n.setSize(mesh_.nCells());
    }

}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simplifiedChapmanEnskog::decompose()
{

    timeSteps_++;

    nAvTimeSteps_++;

    // get cell measurements
    auto& cm = cloud_.cellPropMeasurements();

    forAll(cm.rhoNMean(), iD)
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

    if (timeSteps_ == decompositionInterval_)
    {

        const scalar& deltaT = cloud_.mesh().time().deltaTValue();

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

                if (cloud_.cellCollModel(cell) == cloud_.relCollModel() && cloud_.relaxationCollisionModelName() == "unifiedShakhov")
                {
                    scalar Prandtl = 0.0;
                    scalar viscosity = 0.0;
                    if (translationalT_[cell] > VSMALL)
                    {
                        forAll(typeIds_, iD)
                        {

                            const scalar& mass = cloud_.constProps(iD).mass();
                            const scalar& omega = cloud_.constProps(iD).omega();
                            const scalar& a = cloud_.constProps(iD).alpha();
                            const scalar& d = cloud_.constProps(iD).d();
                            const scalar& rotDoF = cloud_.constProps(iD).rotationalDoF();

                            scalar speciesViscRef = 
                                1.25*(1.0+a)*(2.0+a)*sqrt(mass*physicoChemical::k.value()*cloud_.collTref())
                                /(a*(5.0-2.0*omega)*(7.0-2.0*omega)*sqrt(mathematical::pi)*sqr(d));
                                
                            viscosity += nParcels_[iD][cell]*speciesViscRef*pow(translationalT_[cell]/cloud_.collTref(),omega);

                            Prandtl += nParcels_[iD][cell]*(5.0+rotDoF)/(7.5+rotDoF);

                        }
                        viscosity /= rhoNMean_[cell];
                        Prandtl /= rhoNMean_[cell]; 
                    
                        scalar tau = 0.5*p_[cell]/viscosity*deltaT;
                    
                        heatFluxVector_[cell] /= (1.0 + Prandtl*tau);
                        pressureTensor_[cell] /= (1.0 + tau);
                        shearStressTensor_[cell] /= (1.0 + tau);
                    }
                }

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

        rhoN_.correctBoundaryConditions();
        rhoM_.correctBoundaryConditions();
        p_.correctBoundaryConditions();
        translationalT_.correctBoundaryConditions();
        UMean_.correctBoundaryConditions();
        heatFluxVector_.correctBoundaryConditions();
        pressureTensor_.correctBoundaryConditions();
        shearStressTensor_.correctBoundaryConditions();

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
        scalar instB;
        forAll(mesh_.cells(), cell)
        {

            if (rhoNMean_[cell] > VSMALL && translationalT_[cell] > VSMALL)
            {

                scalar u0(
                    cloud_.maxwellianMostProbableSpeed
                    (
                        translationalT_[cell],
                        rhoMMean_[cell]/rhoNMean_[cell]
                    )
                );

                instB = 0.0;
                forAll(heatFluxVector_[cell], i) 
                {
                    instB = max(instB, fabs(2.0*heatFluxVector_[cell][i]/(p_[cell]*u0)));
                }

                forAll(shearStressTensor_[cell], i) 
                {
                    instB = max(instB, fabs(shearStressTensor_[cell][i]/p_[cell]));
                }
            }
            else
            {
                instB = GREAT;
            }

            // time average macroscopic quantities
            B_[cell] = theta_*instB + (1.0-theta_)*B_[cell];

        }

        B_.correctBoundaryConditions();

        // determine cell collision model
        forAll(mesh_.cells(), cell)
        {
            if (B_[cell] > breakdownMax_)
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

        if (resetAtDecomposition_ && mesh_.time().value() < resetAtDecompositionUntilTime_+0.5*cloud_.mesh().time().deltaTValue())
        {

            nAvTimeSteps_ = 0;

            // reset cell information
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

            forAll(nParcels_, i)
            {
                forAll(nParcels_[i], cell)
                {
                    nParcels_[i][cell] = 0.0;
                }
            }            
        }

        update();

    }

    

}

void Foam::simplifiedChapmanEnskog::update()
{

    // The main properties should be updated first
    updateProperties();

    propsDict_ = hybridDecompositionDict_.subDict(typeName + "Properties");

    propsDict_.readIfPresent("breakdownMax", breakdownMax_);

    propsDict_.readIfPresent("theta", theta_);

    propsDict_.readIfPresent("smoothingPasses", smoothingPasses_);  

}

// ************************************************************************* //

