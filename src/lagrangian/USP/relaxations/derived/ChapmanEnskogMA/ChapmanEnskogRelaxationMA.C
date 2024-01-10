/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "ChapmanEnskogRelaxationMA.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(ChapmanEnskogRelaxationMA, 0);

addToRunTimeSelectionTable(relaxationModel, ChapmanEnskogRelaxationMA, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChapmanEnskogRelaxationMA::ChapmanEnskogRelaxationMA
(
    const dictionary& dict,
    const polyMesh& mesh,
    uspCloud& cloud
)
:
    relaxationModel(dict, mesh, cloud),
    movingAvgCounter_(0),
    movingAvgSpan_(100),
    nTimeSteps_(0),
    infoCounter_(0),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoNMeanInt_(mesh_.nCells(), 0.0),
    molsElec_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    linearKEMeanXnParticle_(mesh_.nCells(), 0.0),
    rotationalEMean_(mesh_.nCells(), 0.0),
    rotationalDofMean_(mesh_.nCells(), 0.0),
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
    momentumMean_(mesh.nCells(), Zero),
    momentumMeanXnParticle_(mesh.nCells(), Zero),
    vibrationalETotal_(),
    electronicETotal_(),
    nParcels_(),
    mccSpecies_(),
    nGroundElectronicLevel_(),
    nFirstElectronicLevel_(),
    boundaryCells_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    rotationalEBF_(),
    rotationalDofBF_(),
    vibTxvDofBF_(),
    totalvDofBF_(),
    speciesRhoNIntBF_(),
    speciesRhoNElecBF_(),
    momentumBF_(),
    vibrationalEBF_(),
    electronicEBF_(),
    speciesRhoNBF_(),
    mccSpeciesBF_(),
    vibTBF_(),
    vDofBF_(),
    Prandtl_
    (
        IOobject
        (
            "Prandtl_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    relaxFreq_
    (
        IOobject
        (
            "relaxFreq_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    ),
    viscosity_
    (
        IOobject
        (
            "viscosity_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure*dimTime, Zero)
    ),
    conductivity_
    (
        IOobject
        (
            "conductivity_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(1, 1, -3, -1, 0), Zero)
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
            "rhoN_",
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
    rotationalT_
    (
        IOobject
        (
            "rotationalT_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, Zero)
    ),
    vibrationalT_
    (
        IOobject
        (
            "vibrationalT_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, Zero)
    ),
    electronicT_
    (
        IOobject
        (
            "electronicT_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, Zero)
    ),
    overallT_
    (
        IOobject
        (
            "overallT_",
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
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimLength/dimTime, Zero)
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
        dimensionedVector(dimensionSet(1, 0, -3, 0, 0), Zero),
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
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero)
    ),
    stressTensor_
    (
        IOobject
        (
            "stressTensor_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    gradT_
    (
        IOobject
        (
            "gradT_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimTemperature/dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    gradU_
    (
        IOobject
        (
            "gradU_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimVelocity/dimLength, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    translationalTAvg_
    (
        IOobject
        (
            "translationalTAvg_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, Zero)
    ),
    UMeanAvg_
    (
        IOobject
        (
            "UMeanAvg_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimVelocity, Zero)
    ),
    heatFluxVectorAvg_
    (
        IOobject
        (
            "heatFluxVectorAvg_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimensionSet(1, 0, -3, 0, 0), Zero)
    ),
    stressTensorAvg_
    (
        IOobject
        (
            "stressTensorAvg_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero)
    ),
    phi_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(sqr(dimLength)*dimVelocity, Zero)
    )
{

    typeIds_.setSize(cloud_.typeIdList().size());
    forAll(typeIds_, iD)
    {
        typeIds_[iD] = iD;
    }

    nGroundElectronicLevel_.setSize(typeIds_.size());

    for (auto& l : nGroundElectronicLevel_)
    {
        l.setSize(mesh_.nCells(), 0.0);
    }

    nFirstElectronicLevel_.setSize(typeIds_.size());

    for (auto& l : nFirstElectronicLevel_)
    {
        l.setSize(mesh_.nCells(), 0.0);
    }

    vibrationalETotal_.setSize(typeIds_.size());

    electronicETotal_.setSize(typeIds_.size());

    for (auto& e : electronicETotal_)
    {
        e.setSize(mesh_.nCells(), 0.0);
    }

    nParcels_.setSize(typeIds_.size());

    for (auto& n : nParcels_)
    {
        n.setSize(mesh_.nCells());
    }

    mccSpecies_.setSize(typeIds_.size());

    for (auto& m : mccSpecies_)
    {
        m.setSize(mesh_.nCells());
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
    rotationalEBF_.setSize(mesh_.boundaryMesh().size());
    rotationalDofBF_.setSize(mesh_.boundaryMesh().size());
    vibTxvDofBF_.setSize(mesh_.boundaryMesh().size());
    totalvDofBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNIntBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNElecBF_.setSize(mesh_.boundaryMesh().size());

    forAll(rhoNBF_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];

        rhoNBF_[j].setSize(patch.size(), 0.0);
        rhoMBF_[j].setSize(patch.size(), 0.0);
        linearKEBF_[j].setSize(patch.size(), 0.0);
        momentumBF_[j].setSize(patch.size(), Zero);
        rotationalEBF_[j].setSize(patch.size(), 0.0);
        rotationalDofBF_[j].setSize(patch.size(), 0.0);
        vibTxvDofBF_[j].setSize(patch.size(), 0.0);
        totalvDofBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNIntBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNElecBF_[j].setSize(patch.size(), 0.0);
    }

    vibrationalEBF_.setSize(typeIds_.size());
    electronicEBF_.setSize(typeIds_.size());
    speciesRhoNBF_.setSize(typeIds_.size());
    mccSpeciesBF_.setSize(typeIds_.size());
    vibTBF_.setSize(typeIds_.size());
    vDofBF_.setSize(typeIds_.size());

    forAll(vibrationalEBF_, i)
    {
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
        speciesRhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());
        vibTBF_[i].setSize(mesh_.boundaryMesh().size());
        vDofBF_[i].setSize(mesh_.boundaryMesh().size());

        forAll(vibrationalEBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
            electronicEBF_[i][j].setSize(patch.size(), 0.0);
            speciesRhoNBF_[i][j].setSize(patch.size(), 0.0);
            mccSpeciesBF_[i][j].setSize(patch.size(), 0.0);
            vibTBF_[i][j].setSize(patch.size(), 0.0);
            vDofBF_[i][j].setSize(patch.size(), 0.0);
        }
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ChapmanEnskogRelaxationMA::calculateProperties()
{

    // Obtain cell measurements

    auto& cm = cloud_.cellPropMeasurements();

    forAll(typeIds_, iD)
    {

        rhoNMean_ += cm.rhoNMean()[iD];
        rhoMMean_ += cm.rhoMMean()[iD];
        linearKEMean_ += cm.linearKEMean()[iD];
        momentumMean_ += cm.momentumMean()[iD];
        rotationalEMean_ += cm.rotationalEMean()[iD];
        rotationalDofMean_ += cm.rotationalDofMean()[iD];
        electronicETotal_[iD] += cm.electronicETotal()[iD];
        rhoNMeanXnParticle_ += cm.rhoNMeanXnParticle()[iD];
        rhoMMeanXnParticle_ += cm.rhoMMeanXnParticle()[iD];
        momentumMeanXnParticle_ += cm.momentumMeanXnParticle()[iD];
        linearKEMeanXnParticle_ += cm.linearKEMeanXnParticle()[iD];

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

        rhoNMeanInt_ += cm.rhoNMeanInt()[iD];
        molsElec_ += cm.molsElec()[iD];

        nParcels_[iD] += cm.nParcels()[iD];
        mccSpecies_[iD] += cm.mccSpecies()[iD];

        nGroundElectronicLevel_[iD] += cm.nGroundElectronicLevel()[iD];
        nFirstElectronicLevel_[iD] += cm.nFirstElectronicLevel()[iD];

        forAll(vibrationalETotal_[iD], v)
        {
            vibrationalETotal_[iD][v] += cm.vibrationalETotal()[iD][v];
        }

    }

    // Obtain boundary measurements
    const scalar nParticle = cloud_.nParticle();
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
                rotationalEBF_[j][k] += bm.rotationalEBF()[i][j][k];
                rotationalDofBF_[j][k] += bm.rotationalDofBF()[i][j][k];
                speciesRhoNBF_[iD][j][k] += bm.rhoNBF()[i][j][k];
                vibrationalEBF_[iD][j][k] += bm.vibrationalEBF()[i][j][k];
                electronicEBF_[iD][j][k] += bm.electronicEBF()[i][j][k];
                mccSpeciesBF_[iD][j][k] += bm.mccSpeciesBF()[i][j][k];
                speciesRhoNIntBF_[j][k] += bm.rhoNIntBF()[i][j][k];
                speciesRhoNElecBF_[j][k] += bm.rhoNElecBF()[i][j][k];
            }
        }
        
    }

    // Calculate internal vol fields
    const scalar nAvTimeSteps = scalar(nTimeSteps_);

    forAll(mesh_.cells(), cell)
    {

        if (rhoNMean_[cell] > VSMALL)
        {

            const scalar cellVolume = mesh_.cellVolumes()[cell];

            // Density
            rhoN_[cell] =
                (rhoNMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
            rhoM_[cell] =
                (rhoMMeanXnParticle_[cell])/(nAvTimeSteps*cellVolume);
            
            // Velcoity
            scalar rhoMMean =
                rhoMMeanXnParticle_[cell]/(cellVolume*nAvTimeSteps);
            UMean_[cell] =
                momentumMeanXnParticle_[cell]
               /(rhoMMean*cellVolume*nAvTimeSteps);
            scalar linearKEMean =
                0.5*linearKEMeanXnParticle_[cell]
               /(cellVolume*nAvTimeSteps);
            scalar rhoNMean =
                rhoNMeanXnParticle_[cell]/(cellVolume*nAvTimeSteps);

            // Translational temperature
            translationalT_[cell] =
                2.0/(3.0*physicoChemical::k.value()*rhoNMean)
               *(
                    linearKEMean
                  - 0.5*rhoMMean*(UMean_[cell] & UMean_[cell])
                );

            // Pressure
            p_[cell] =
                rhoN_[cell]*physicoChemical::k.value()
                *translationalT_[cell];

            // Pressure tensor
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

            scalar scalarPressure = (1.0/3.0)*
                                    (pressureTensor_[cell].xx() +
                                    pressureTensor_[cell].yy() +
                                    pressureTensor_[cell].zz());

            // Stress tensor
            stressTensor_[cell] = -pressureTensor_[cell];
            stressTensor_[cell].xx() += scalarPressure;
            stressTensor_[cell].yy() += scalarPressure;
            stressTensor_[cell].zz() += scalarPressure;

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

            // Relaxation frequency !!!Check mixtures and vibrational-electronic DoF
            Prandtl_[cell] = 0.0;
            viscosity_[cell] = 0.0;
            conductivity_[cell] = 0.0;
            List<scalar> speciesPrandtl(typeIds_.size(), Zero);
            List<scalar> speciesVisc(typeIds_.size(), Zero);
            List<scalar> speciesCond(typeIds_.size(), Zero);

            if (translationalT_[cell] > VSMALL)
            {

                forAll(typeIds_, iD)
                {

                    const scalar& Tref = cloud_.constProps(iD).Tref();
                    const scalar& mass = cloud_.constProps(iD).mass();
                    const scalar& omega = cloud_.constProps(iD).omega();
                    const scalar& d = cloud_.constProps(iD).d();
                    const scalar& rotDoF = cloud_.constProps(iD).rotationalDoF();


                    speciesPrandtl[iD] = (5+rotDoF)/(7.5+rotDoF);
                    Prandtl_[cell] += nParcels_[iD][cell]*speciesPrandtl[iD];

                    scalar speciesViscRef = 
                        7.5*sqrt(mass*physicoChemical::k.value()*Tref)
                        /(sqrt(mathematical::pi)*(5.0-2.0*omega)*(7.0-2.0*omega)*sqr(d));
                    speciesVisc[iD] = speciesViscRef*pow(translationalT_[cell]/Tref,omega);
                    viscosity_[cell] += nParcels_[iD][cell]*speciesVisc[iD];

                    speciesCond[iD] = ((2.5+0.5*rotDoF)/speciesPrandtl[iD])*(physicoChemical::k.value()/mass)*speciesVisc[iD];
                    conductivity_[cell] += nParcels_[iD][cell]*speciesCond[iD];
 
                }
                viscosity_[cell] /= rhoNMean_[cell];
                Prandtl_[cell] /= rhoNMean_[cell];
                conductivity_[cell] /= rhoNMean_[cell];

                relaxFreq_[cell] = p_[cell]/viscosity_[cell];

            }
            else
            {
                relaxFreq_[cell] = 0.0;
            }

        }
        else
        {
            viscosity_[cell] = 0.0;
            conductivity_[cell] = 0.0;
            rhoN_[cell] = 0.0;
            p_[cell] = 0.0;
            UMean_[cell] = vector::zero;
            translationalT_[cell] = 0.0;
            pressureTensor_[cell] = tensor::zero;
            stressTensor_[cell] = tensor::zero;
            heatFluxVector_[cell] =vector::zero;
        }

    }

    // Calcualte boundary vol fields
    List<scalarField> vibTBF(mesh_.boundaryMesh().size());

    forAll(rhoNBF_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];

        vibTBF[j].setSize(patch.size(), 0.0);

        if (isA<wallPolyPatch>(patch))
        {
            forAll(rhoN_.boundaryFieldRef()[j], k)
            {
                rhoN_.boundaryFieldRef()[j][k] = rhoNBF_[j][k]*nParticle/nAvTimeSteps;
                rhoM_.boundaryFieldRef()[j][k] = rhoMBF_[j][k]*nParticle/nAvTimeSteps;

                if (rhoM_.boundaryFieldRef()[j][k] > VSMALL)
                {
                    UMean_.boundaryFieldRef()[j][k] =
                        momentumBF_[j][k]*nParticle
                       /(rhoM_.boundaryFieldRef()[j][k]*nAvTimeSteps);
                }
                else
                {
                    UMean_.boundaryFieldRef()[j][k] = vector::zero;
                }

                scalar rhoMMean = rhoMBF_[j][k]*nParticle/nAvTimeSteps;
                scalar linearKEMean = linearKEBF_[j][k]*nParticle/nAvTimeSteps;
                scalar rhoNMean = rhoNBF_[j][k]*nParticle/nAvTimeSteps;

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

                if (rotationalDofBF_[j][k] > VSMALL)
                {
                    rotationalT_.boundaryFieldRef()[j][k] =
                        (2.0/physicoChemical::k.value())
                       *(rotationalEBF_[j][k]/rotationalDofBF_[j][k]);
                }
                else
                {
                    rotationalT_.boundaryFieldRef()[j][k] = 0.0;
                }

                // electronic temperature
                scalar totalEDof = 0.0;
                scalar elecT = 0.0;

                electronicT_.boundaryFieldRef()[j][k] = elecT;

                scalar nRotDof = 0.0;

                if (rhoNBF_[j][k] > VSMALL)
                {
                    nRotDof = rotationalDofBF_[j][k]/rhoNBF_[j][k];
                }

                overallT_.boundaryFieldRef()[j][k] =
                    (
                        (3.0*translationalT_.boundaryFieldRef()[j][k])
                      + (nRotDof*rotationalT_.boundaryFieldRef()[j][k])
                      + (
                            totalvDofBF_[j][k]
                           *vibrationalT_.boundaryFieldRef()[j][k]
                        )
                      + (totalEDof*elecT)
                    )
                   /(3.0 + nRotDof + totalvDofBF_[j][k] + totalEDof);

                totalvDofBF_[j][k] = 0.0;

            }

        }
    }


    translationalT_[0] = 301.8116398017609;
    translationalT_[mesh_.nCells()] = 301.8116398017609;

    UMean_[0] = vector(195.6925272210719, 0, 0);
    UMean_[mesh_.nCells()-1] = vector(195.6925272210719, 0, 0);

    forAll(mesh_.boundary(), patchI)
    {
    
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];
    
        // Loop over all faces of boundary patch
        forAll(mesh_.boundary()[patchI], faceI)
        {
            const label& bCell = mesh_.boundary()[patchI].faceCells()[faceI];    // Boundary cell index
            const label& bFace = mesh_.boundary()[patchI].start() + faceI;        // Face index
                
            if
            (
                isA<polyPatch>(patch)
            && !isA<emptyPolyPatch>(patch)
            && !isA<cyclicPolyPatch>(patch)
            )
            {
    
                if (patchI==2)
                {
                    translationalT_.boundaryFieldRef()[patchI][faceI] = 301.5572472698444;
                    UMean_.boundaryFieldRef()[patchI][faceI] = vector(196.7025409516354, 0, 0);
                }
                else
                {
                    translationalT_.boundaryFieldRef()[patchI][faceI] = 301.5572472698444;
                    UMean_.boundaryFieldRef()[patchI][faceI] = vector(-196.7025409516354, 0, 0);
                }
                
    
            }
        }
    
    }

    /*forAll(mesh_.boundary(), patchI)
    {

        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        // Loop over all faces of boundary patch
        forAll(mesh_.boundary()[patchI], faceI)
        {
            const label& bCell = mesh_.boundary()[patchI].faceCells()[faceI];    // Boundary cell index
            const label& bFace = mesh_.boundary()[patchI].start() + faceI;        // Face index
                
            if
            (
                isA<polyPatch>(patch)
            && !isA<emptyPolyPatch>(patch)
            && !isA<cyclicPolyPatch>(patch)
            )
            {

                gradT_[bCell] = vector::zero;
                gradU_[bCell] = tensor::zero;
                const vector& cellCentre = mesh_.cellCentres()[bCell];
                const vector& faceCentre = mesh_.faceCentres()[bFace];
                forAll(mesh_.cellCells()[bCell], i)
                {
                    const label& adjCell = mesh_.cellCells()[bCell][i];
                    const vector directionVector = mesh_.cellCentres()[bCell]-mesh_.cellCentres()[adjCell];
                    const scalar directionVectorMag = mag(directionVector);
                    gradT_[bCell] += directionVector/sqr(directionVectorMag)*(translationalT_[bCell] - translationalT_[adjCell]);
                    gradU_[bCell] += directionVector/sqr(directionVectorMag)*(UMean_[bCell] - UMean_[adjCell]);
                }
                gradT_[bCell] /= mesh_.cellCells()[bCell].size();
                gradU_[bCell] /= mesh_.cellCells()[bCell].size();

                translationalT_.boundaryFieldRef()[patchI][faceI] = translationalT_[bCell] + (gradT_[bCell] & (faceCentre - cellCentre));
                UMean_.boundaryFieldRef()[patchI][faceI] = UMean_[bCell] + (gradU_[bCell] & (faceCentre - cellCentre));

            }
        }

    }*/

    /*forAll(mesh_.boundary(), patchI)
    {

        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        // Loop over all faces of boundary patch
        forAll(mesh_.boundary()[patchI], faceI)
        {
            const label& bCell = mesh_.boundary()[patchI].faceCells()[faceI];    // Boundary cell index
            const label& bFace = mesh_.boundary()[patchI].start() + faceI;        // Face index
                
            if
            (
                isA<polyPatch>(patch)
            && !isA<emptyPolyPatch>(patch)
            && !isA<cyclicPolyPatch>(patch)
            )
            {

                translationalT_.boundaryFieldRef()[patchI][faceI] = translationalT_[bCell];
                UMean_.boundaryFieldRef()[patchI][faceI] = UMean_[bCell];

            }
        }

    }*/

    /*forAll(mesh_.boundary(), patchI)
    {

        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        // Loop over all faces of boundary patch
        forAll(mesh_.boundary()[patchI], faceI)
        {
            const label& bCell = mesh_.boundary()[patchI].faceCells()[faceI];    // Boundary cell index
            const label& bFace = mesh_.boundary()[patchI].start() + faceI;        // Face index
                
            if
            (
                isA<polyPatch>(patch)
            && !isA<emptyPolyPatch>(patch)
            && !isA<cyclicPolyPatch>(patch)
            )
            {

                gradT_[bCell] = vector::zero;
                gradU_[bCell] = tensor::zero;
                const vector& cellCentre = mesh_.cellCentres()[bCell];
                const vector& faceCentre = mesh_.faceCentres()[bFace];
                forAll(mesh_.cellCells()[bCell], i)
                {
                    const label& adjCell = mesh_.cellCells()[bCell][i];
                    const vector directionVector = mesh_.cellCentres()[bCell]-mesh_.cellCentres()[adjCell];
                    const scalar directionVectorMag = mag(directionVector);
                    gradT_[bCell] += directionVector/sqr(directionVectorMag)*(translationalT_[bCell] - translationalT_[adjCell]);
                    gradU_[bCell] += directionVector/sqr(directionVectorMag)*(UMean_[bCell] - UMean_[adjCell]);
                }
                gradT_[bCell] /= mesh_.cellCells()[bCell].size();
                gradU_[bCell] /= mesh_.cellCells()[bCell].size();

                translationalT_.boundaryFieldRef()[patchI][faceI] = translationalT_[bCell] + (gradT_[bCell] & (faceCentre - cellCentre));
                UMean_.boundaryFieldRef()[patchI][faceI] = UMean_[bCell] + (gradU_[bCell] & (faceCentre - cellCentre));

            }
        }

    }

    //phi_ = fvc::flux(UMean_);

    forAll(mesh_.boundary(), patchI)
    {

        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        // Loop over all faces of boundary patch
        forAll(mesh_.boundary()[patchI], faceI)
        {

            const label& bCell = mesh_.boundary()[patchI].faceCells()[faceI];    // Boundary cell index

            gradT_[bCell] = vector::zero;
            gradU_[bCell] = tensor::zero;
            forAll(mesh_.cellCells()[bCell], i)
            {
                const label& adjCell = mesh_.cellCells()[bCell][i];
                const vector directionVector = mesh_.cellCentres()[bCell]-mesh_.cellCentres()[adjCell];
                const scalar directionVectorMag = mag(directionVector);
                gradT_[bCell] += directionVector/sqr(directionVectorMag)*(translationalT_[bCell] - translationalT_[adjCell]);
                gradU_[bCell] += directionVector/sqr(directionVectorMag)*(UMean_[bCell] - UMean_[adjCell]);
            }
            gradT_[bCell] /= mesh_.cellCells()[bCell].size();
            gradU_[bCell] /= mesh_.cellCells()[bCell].size();

        }

    }*/

    //std::cout << translationalT_.boundaryField()[3][0] << " " << translationalT_[0] << " " << translationalT_[1] << " " << std::endl;

    //movingAvgCounter_++;
    //translationalTAvg_ = 
    //    ((min(movingAvgCounter_,movingAvgSpan_)-1.0)*translationalTAvg_ + translationalT_)/min(movingAvgCounter_,movingAvgSpan_);
    //UMeanAvg_ = 
    //    ((min(movingAvgCounter_,movingAvgSpan_)-1.0)*UMeanAvg_ + UMean_)/min(movingAvgCounter_,movingAvgSpan_);

    //gradT_ = fvc::grad(translationalTAvg_);
    //gradU_ = fvc::grad(UMeanAvg_);


    /*forAll(mesh_.cells(), cell)
    {

        // initialize neighborCells
        DynamicList<label> neighborCells;
        neighborCells.setCapacity(0);

        // get cells in extended neighbor area
        label neighborCandidate;
        label neighborCellinitial = 0;
        label neighborCellFinal = 0;
        label neighborCellCapacity = 1;
        boolList isNeighbor(mesh_.nCells());
        isNeighbor = false;
        isNeighbor[cell]=true;
        neighborCells.append(cell);

        for(label level=1; level<=2; level++)
        {

            neighborCellinitial=neighborCellFinal;
            neighborCellFinal=neighborCellCapacity-1;

            for(label i=neighborCellinitial; i<=neighborCellFinal; i++)
            {

                forAll(mesh_.cellCells()[neighborCells[i]], j)
                {

                    neighborCandidate=mesh_.cellCells()[neighborCells[i]][j];

                    if (!isNeighbor[neighborCandidate])
                    {

                        neighborCellCapacity++;
                        isNeighbor[neighborCandidate]=true;
                        neighborCells.append(neighborCandidate);

                    }
                }
            }
        }

        //if (cell == 0 || cell == mesh_.nCells()-1)
        //{
        //    std::cout << cell << " " << neighborCells[1] << " " << neighborCells[2] << std::endl;
        //}
        //else
        //{
        //    std::cout << cell << " " << neighborCells[1] << " " << neighborCells[2] << " " << neighborCells[3] << " " << neighborCells[4] << std::endl;
        //}

        gradT_[cell] = vector::zero;
        gradU_[cell] = tensor::zero;
        const vector& cellCentre = mesh_.cellCentres()[cell];
        scalar directionVectorMagSum = 0.0;
        for(label i = 1; i<=neighborCellFinal; i++)
        {
            const label& nCell = neighborCells[i];
            const vector directionVector = mesh_.cellCentres()[cell]-mesh_.cellCentres()[nCell];
            const scalar directionVectorMag = mag(directionVector);
            directionVectorMagSum += directionVectorMag;
            gradT_[cell] += directionVector/directionVectorMag*(translationalT_[cell] - translationalT_[nCell]);
            gradU_[cell] += directionVector/directionVectorMag*(UMean_[cell] - UMean_[nCell]);
        }
        gradT_[cell] /= directionVectorMagSum;
        gradU_[cell] /= directionVectorMagSum;
        //std::cout << cell << " " << gradT_[cell].y() << std::endl;
    }*/

    /*forAll(mesh_.cells(), cell)
    {

        // initialize neighborCells
        DynamicList<label> neighborCells;
        neighborCells.setCapacity(0);

        // get cells in extended neighbor area
        label neighborCandidate;
        label neighborCellinitial = 0;
        label neighborCellFinal = 0;
        label neighborCellCapacity = 1;
        boolList isNeighbor(mesh_.nCells());
        isNeighbor = false;
        isNeighbor[cell]=true;
        neighborCells.append(cell);

        for(label level=1; level<=2; level++)
        {

            neighborCellinitial=neighborCellFinal;
            neighborCellFinal=neighborCellCapacity-1;

            for(label i=neighborCellinitial; i<=neighborCellFinal; i++)
            {

                forAll(mesh_.cellCells()[neighborCells[i]], j)
                {

                    neighborCandidate=mesh_.cellCells()[neighborCells[i]][j];

                    if (!isNeighbor[neighborCandidate])
                    {

                        neighborCellCapacity++;
                        isNeighbor[neighborCandidate]=true;
                        neighborCells.append(neighborCandidate);

                    }
                }
            }
        }

        //if (cell == 0 || cell == mesh_.nCells()-1)
        //{
        //    std::cout << cell << " " << neighborCells[1] << " " << neighborCells[2] << std::endl;
        //}
        //else
        //{
        //    std::cout << cell << " " << neighborCells[1] << " " << neighborCells[2] << " " << neighborCells[3] << " " << neighborCells[4] << std::endl;
        //}

        gradT_[cell] = vector::zero;
        gradU_[cell] = tensor::zero;
        const vector& cellCentre = mesh_.cellCentres()[cell];
        scalar directionVectorMagSum = 0.0;
        for(label i = 1; i<=neighborCellFinal; i++)
        {
            const label& nCell = neighborCells[i];
            const vector directionVector = mesh_.cellCentres()[cell]-mesh_.cellCentres()[nCell];
            const scalar directionVectorMag = mag(directionVector);
            directionVectorMagSum += directionVectorMag;
            gradT_[cell] += directionVector/directionVectorMag*(translationalT_[cell] - translationalT_[nCell]);
            gradU_[cell] += directionVector/directionVectorMag*(UMean_[cell] - UMean_[nCell]);
        }
        gradT_[cell] /= directionVectorMagSum;
        gradU_[cell] /= directionVectorMagSum;
        //std::cout << cell << " " << gradT_[cell].y() << std::endl;
    }*/


    /*forAll(mesh_.cells(), cellI)
    {

        gradT_[cellI] = vector::zero;
        gradU_[cellI] = tensor::zero;
        forAll(mesh_.cellCells()[cellI], i)
        {
            const label& cellJ = mesh_.cellCells()[cellI][i];
            const vector directionVector = mesh_.cellCentres()[cellI]-mesh_.cellCentres()[cellJ];
            const scalar directionVectorMag = mag(directionVector);
            gradT_[cellI] += directionVector/sqr(directionVectorMag)*(translationalT_[cellI] - translationalT_[cellJ]);
            gradU_[cellI] += directionVector/sqr(directionVectorMag)*(UMean_[cellI] - UMean_[cellJ]);
        }
        gradT_[cellI] /= mesh_.cellCells()[cellI].size();
        gradU_[cellI] /= mesh_.cellCells()[cellI].size();
        //std::cout << cell << " " << gradT_[cell].y() << std::endl;
    }*/

}

void Foam::ChapmanEnskogRelaxationMA::resetProperties()
{

    if (nTimeSteps_ == 1)
    {
        nTimeSteps_ = 0;

        forAll(mesh_.cells(), cell)
        {
            rhoNMean_[cell] = scalar(0.0);
            rhoMMean_[cell] = scalar(0.0);
            linearKEMean_[cell] = scalar(0.0);
            momentumMean_[cell] = vector::zero;
            rotationalEMean_[cell] = scalar(0.0);
            rotationalDofMean_[cell] = scalar(0.0);
            rhoNMeanInt_[cell] = scalar(0.0);
            molsElec_[cell] = scalar(0.0),
            muu_[cell] = scalar(0.0);
            muv_[cell] = scalar(0.0);
            muw_[cell] = scalar(0.0);
            mvv_[cell] = scalar(0.0);
            mvw_[cell] = scalar(0.0);
            mww_[cell] = scalar(0.0);
            mcc_[cell] = scalar(0.0);
            mccu_[cell] = scalar(0.0);
            mccv_[cell] = scalar(0.0);
            mccw_[cell] = scalar(0.0);
            eu_[cell] = scalar(0.0);
            ev_[cell] = scalar(0.0);
            ew_[cell] = scalar(0.0);
            e_[cell] = scalar(0.0);
            rhoNMeanXnParticle_[cell] = scalar(0.0);
            rhoMMeanXnParticle_[cell] = scalar(0.0);
            momentumMeanXnParticle_[cell] = vector::zero;
            linearKEMeanXnParticle_[cell] = scalar(0.0);

            forAll(electronicETotal_, iD)
            {

                electronicETotal_[iD][cell] = 0.0;
                mccSpecies_[iD][cell] = 0.0;
                nParcels_[iD][cell] = 0.0;
                nGroundElectronicLevel_[iD][cell] = 0.0;
                nFirstElectronicLevel_[iD][cell] = 0.0;

                forAll(vibrationalETotal_[iD], v)
                {
                    vibrationalETotal_[iD][v][cell] = 0.0;
                }
            }
        }

         // reset boundary information
        forAll(rhoNBF_, j)
        {
            rhoNBF_[j] = 0.0;
            rhoMBF_[j] = 0.0;
            linearKEBF_[j] = 0.0;
            speciesRhoNIntBF_[j] = 0.0;
            speciesRhoNElecBF_[j] = 0.0;
            rotationalEBF_[j] = 0.0;
            rotationalDofBF_[j] = 0.0;
            momentumBF_[j] = vector::zero;
        }

        forAll(speciesRhoNBF_, i)
        {
            forAll(speciesRhoNBF_[i], j)
            {
                speciesRhoNBF_[i][j] = 0.0;
                vibrationalEBF_[i][j] = 0.0;
                electronicEBF_[i][j] = 0.0;
                mccSpeciesBF_[i][j] = 0.0;
            }
        }

    }

}

void Foam::ChapmanEnskogRelaxationMA::relax()
{

    nTimeSteps_++;

    const scalar& deltaT = cloud_.mesh().time().deltaTValue();
    const scalar& nParticle = cloud_.nParticle();
    
    // Relax particles to local distribution
    label relaxations = 0;

    List<DynamicList<uspParcel*>>&
        cellOccupancy = cloud_.cellOccupancy();

    // Calculate required macroscopic properties
    calculateProperties();

            //pressureTensor_[cell].xx() = rhoN_[cell]*
            //(
            //    muu_[cell]/(rhoNMean_[cell]) -
            //    (
            //        (rhoMMean_[cell]/(rhoNMean_[cell]))
            //        *UMean_[cell].x()*UMean_[cell].x()
            //    )
            //);

    //std::cout << "CE " << rhoN_[0] << " " << muv_[0] << " " << rhoNMean_[0] << " " << rhoMMean_[0] << " " << UMean_[0].x() << " " << UMean_[0].y() << " " << stressTensor_[0].xy() << std::endl;

    //std::ofstream outfile;
//
    //outfile.open("./T.dat",std::ios_base::app);
    //forAll(cellOccupancy, cell)
    //{
    //    outfile << translationalT_[cell] << " ";
    //}
    //outfile << std::endl;
    //outfile.close();
//
    //outfile.open("./p.dat",std::ios_base::app);
    //forAll(cellOccupancy, cell)
    //{
    //    outfile << p_[cell] << " ";
    //}
    //outfile << std::endl;
    //outfile.close();
//
    //outfile.open("./Pxy.dat",std::ios_base::app);
    //forAll(cellOccupancy, cell)
    //{
    //    outfile << stressTensor_[cell].xy() << " ";
    //}
    //outfile << std::endl;
    //outfile.close();
//
    //outfile.open("./Qy.dat",std::ios_base::app);
    //forAll(cellOccupancy, cell)
    //{
    //    outfile << heatFluxVector_[cell].y() << " ";
    //}
    //outfile << std::endl;
    //outfile.close();

    /*translationalTAvg_ = translationalT_;
    UMeanAvg_ = UMean_;
    for (label pass=0; pass<10; ++pass)
    {
        translationalTAvg_ = fvc::average(fvc::interpolate(translationalTAvg_));
        UMeanAvg_ = fvc::average(fvc::interpolate(UMeanAvg_));
    }*/

    phi_ = fvc::flux(UMean_);

    gradT_ = fvc::grad(translationalT_);
    gradU_ = fvc::grad(UMean_);

    heatFluxVector_ = -conductivity_*gradT_;
    stressTensor_ = viscosity_*(gradU_ + gradU_.T()-(2e0/3e0)*I*tr(gradU_));

    //for (label pass=0; pass<10; ++pass)
    //{
    //    heatFluxVector_ = fvc::average(fvc::interpolate(heatFluxVector_));
    //    stressTensor_ = fvc::average(fvc::interpolate(stressTensor_));
    //    heatFluxVector_.correctBoundaryConditions();
    //    stressTensor_.correctBoundaryConditions();
    //}

    //outfile.open("./PxyGrad.dat",std::ios_base::app);
    //forAll(cellOccupancy, cell)
    //{
    //    outfile << stressTensor_[cell].xy() << " ";
    //}
    //outfile << std::endl;
    //outfile.close();
//
    //outfile.open("./QyGrad.dat",std::ios_base::app);
    //forAll(cellOccupancy, cell)
    //{
    //    outfile << heatFluxVector_[cell].y() << " ";
    //}
    //outfile << std::endl;
    //outfile.close();

    //for (label pass=0; pass<10; ++pass)
    //{
    //    heatFluxVector_ = fvc::average(fvc::interpolate(heatFluxVector_));
    //    stressTensor_ = fvc::average(fvc::interpolate(stressTensor_));
    //    heatFluxVector_.correctBoundaryConditions();
    //    stressTensor_.correctBoundaryConditions();
    //}

    //Calculate heatflux and stress
    /*movingAvgCounter_++;
    translationalTAvg_ = 
        ((min(movingAvgCounter_,movingAvgSpan_)-1.0)*translationalTAvg_ + translationalT_)/min(movingAvgCounter_,movingAvgSpan_);
    UMeanAvg_ = 
        ((min(movingAvgCounter_,movingAvgSpan_)-1.0)*UMeanAvg_ + UMean_)/min(movingAvgCounter_,movingAvgSpan_);

    heatFluxVectorAvg_ = 
        ((min(movingAvgCounter_,movingAvgSpan_)-1.0)*heatFluxVectorAvg_ + heatFluxVector_)/min(movingAvgCounter_,movingAvgSpan_);
    stressTensorAvg_ = 
        ((min(movingAvgCounter_,movingAvgSpan_)-1.0)*stressTensorAvg_ + stressTensor_)/min(movingAvgCounter_,movingAvgSpan_);

    gradT_ = fvc::grad(translationalTAvg_);
    heatFluxVector_ = -conductivity_*gradT_;
    gradU_ = fvc::grad(UMeanAvg_);
    stressTensor_ = viscosity_*(gradU_ + gradU_.T()-(2e0/3e0)*I*tr(gradU_));

    for (label pass=0; pass<10; ++pass)
    {
        heatFluxVector_ = fvc::average(fvc::interpolate(heatFluxVector_));
        stressTensor_ = fvc::average(fvc::interpolate(stressTensor_));
        heatFluxVector_.correctBoundaryConditions();
        stressTensor_.correctBoundaryConditions();
    }*/

    forAll(cellOccupancy, cell)
    {

        if (cloud_.cellCollModel(cell) == cloud_.relCollModel() && cellOccupancy[cell].size() >= 2)
        {

            const DynamicList<uspParcel*>& cellParcels(cellOccupancy[cell]);

            forAll(cellParcels, i)
            {

                //if (cloud_.rndGen().sample01<scalar>() < 1.0-exp(-relaxFreq_[cell]*deltaT)) //cloud_.rndGen().sample01<scalar>()
                //{

                    uspParcel& parcel = *cellParcels[i];
                    const scalar mass = cloud_.constProps(parcel.typeId()).mass();

                    scalar u0(
                        cloud_.maxwellianMostProbableSpeed
                        (
                            translationalT_[cell],
                            mass
                        )
                    );

                    scalar maxHeatFlux = VSMALL;
                    forAll(heatFluxVector_[cell],i) 
                    {
                        if (maxHeatFlux < fabs(heatFluxVector_[cell][i])) 
                        {
                            maxHeatFlux = fabs(heatFluxVector_[cell][i]);
                        }
                    }
                    maxHeatFlux = 2.0*maxHeatFlux/(p_[cell]*u0);

                    scalar maxStress = VSMALL;
                    forAll(stressTensor_[cell],i) 
                    {
                        if (maxStress < fabs(stressTensor_[cell][i])) 
                        {
                            maxStress = fabs(stressTensor_[cell][i]);
                        }
                    }
                    maxStress = maxStress/p_[cell];

                    scalar breakdownParameter = max(maxHeatFlux,maxStress);

                    // Relax particle
                    parcel.U() = samplePostRelaxationVelocity
                            (
                                breakdownParameter,
                                mass,
                                p_[cell],
                                translationalT_[cell],
                                UMean_[cell],
                                heatFluxVector_[cell],
                                stressTensor_[cell]
                            );

                    // Relax particle
                    /*parcel.U() = samplePostRelaxationVelocityTest
                            (
                                mass,
                                Prandtl_[cell],
                                p_[cell],
                                translationalT_[cell],
                                UMean_[cell],
                                pressureTensor_[cell]
                            );*/

                    relaxations++;

                //}

            }

            conserveMomentumAndEnergy(cell);

        }
    }    

    resetProperties();

    infoCounter_++;
    reduce(relaxations, sumOp<label>());

    if (infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (relaxations>0)
        {
            Info<< "    Relaxations                      = "
                << relaxations << nl
                << endl;
            infoCounter_ = 0;
        }
        else
        {
            Info<< "    No relaxations" << endl;
            infoCounter_ = 0;
        }
    }

}

void Foam::ChapmanEnskogRelaxationMA::conserveMomentumAndEnergy
(
    const label& cell
)
{

    // Momentum and energy conservation scheme
    const scalar nParticle = cloud_.nParticle();
    const DynamicList<uspParcel*>& cellParcels(cloud_.cellOccupancy()[cell]);

    scalar linearKEMeanXnParticle = 0.0;
    vector momentumMeanXnParticle = vector::zero;

    scalar rhoNMeanXnParticle = 0.0;
    scalar rhoMMeanXnParticle = 0.0;
    scalar postTranslationalT = 0.0;
    vector postUMean = vector::zero;

    forAll(cellParcels, i)
    {

        const uspParcel& p = *cellParcels[i];

        const scalar mass = cloud_.constProps(p.typeId()).mass();
        const vector& U = p.U();
        const scalar& CWF = cloud_.cellWF(p.cell());
        const scalar& RWF = cloud_.axiRWF(p.position());

        rhoNMeanXnParticle += CWF*RWF*nParticle;
        rhoMMeanXnParticle += mass*CWF*RWF*nParticle;
        linearKEMeanXnParticle += mass*(U & U)*CWF*RWF*nParticle;
        momentumMeanXnParticle += mass*U*CWF*RWF*nParticle;

    }

    // Velocity
    postUMean = momentumMeanXnParticle/rhoMMeanXnParticle;

    // Translational temperature
    postTranslationalT =
        1.0/(3.0*physicoChemical::k.value()*rhoNMeanXnParticle)
       *(
            linearKEMeanXnParticle
          - rhoMMeanXnParticle*(postUMean & postUMean)
        );

    forAll(cellParcels, i)
    {
        uspParcel& p = *cellParcels[i];
        p.U() = UMean_[cell] + (p.U()-postUMean)*sqrt(translationalT_[cell]/postTranslationalT);
    }

}

Foam::vector Foam::ChapmanEnskogRelaxationMA::samplePostRelaxationVelocity
(   
    const scalar& B,
    const scalar& m,
    const scalar& p,
    const scalar& T,
    const vector& U,
    const vector& q,
    const tensor& s
)
{

    vector v;
    scalar gamma;
    scalar A = 1.0 + 30.0*B;
    scalar u0(cloud_.maxwellianMostProbableSpeed(T,m));

    do {

        v = cloud_.rndGen().GaussNormal<vector>()/sqrt(2.0);
        gamma = 1.0+2.0/(p*u0)*(q.x()*v.x()+q.y()*v.y()+q.z()*v.z())*((sqr(v.x())+sqr(v.y())+sqr(v.z()))/2.5-1.0)
               -1.0/p*(s.xx()*sqr(v.x())+s.yy()*sqr(v.y())+s.zz()*sqr(v.z())+2.0*s.xy()*v.x()*v.y()+2.0*s.xz()*v.x()*v.z()+2.0*s.yz()*v.y()*v.z());

    } while(A*cloud_.rndGen().sample01<scalar>()>gamma);

    return U+u0*v;

}

Foam::vector Foam::ChapmanEnskogRelaxationMA::samplePostRelaxationVelocityTest
(   
    const scalar& m,
    const scalar& Pr,
    const scalar& p,
    const scalar& T,
    const vector& U,
    const tensor& pT
)
{

// Sample particle velocity from Maxwellian distribution
vector uP = cloud_.rndGen().GaussNormal<vector>()/sqrt(2.0);

// Calculate approximate transformation tensor for ESBGK model
tensor S = I-0.5*(1-Pr)/Pr*(pT/p-I);

return cloud_.maxwellianMostProbableSpeed(T,m)*(S & uP) + U;

}

// ************************************************************************* //

