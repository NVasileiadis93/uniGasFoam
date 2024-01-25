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

#include "ShakhovRelaxation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(ShakhovRelaxation, 0);

addToRunTimeSelectionTable(relaxationModel, ShakhovRelaxation, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ShakhovRelaxation::ShakhovRelaxation
(
    const dictionary& dict,
    const polyMesh& mesh,
    uspCloud& cloud
)
:
    relaxationModel(dict, mesh, cloud),
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
    momentumMean_(mesh_.nCells(), vector::zero),
    momentumMeanXnParticle_(mesh_.nCells(), vector::zero),
    vibrationalETotal_(),
    electronicETotal_(),
    nParcels_(),
    nParcelsXnParticle_(),
    mccSpecies_(),
    nGroundElectronicLevel_(),
    nFirstElectronicLevel_(),
    Prandtl_
    (
        IOobject
        (
            "Prandtl_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    relaxFreq_
    (
        IOobject
        (
            "relaxFreq_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero),
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
        dimensionedScalar(dimTemperature, Zero),
        zeroGradientFvPatchScalarField::typeName
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
        dimensionedScalar(dimTemperature, Zero),
        zeroGradientFvPatchScalarField::typeName
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
        dimensionedScalar(dimTemperature, Zero),
        zeroGradientFvPatchScalarField::typeName
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
        dimensionedVector(dimVelocity, vector::zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    heatFluxVector_
    (
        IOobject
        (
            "heatFluxVector_",
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
            "pressureTensor_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, tensor::zero),
        zeroGradientFvPatchScalarField::typeName
    )
{

    typeIds_.setSize(cloud_.typeIdList().size());
    forAll(typeIds_, iD)
    {
        typeIds_[iD] = iD;
    }

    nGroundElectronicLevel_.setSize(typeIds_.size());

    nFirstElectronicLevel_.setSize(typeIds_.size());

    electronicETotal_.setSize(typeIds_.size());

    nParcels_.setSize(typeIds_.size());

    nParcelsXnParticle_.setSize(typeIds_.size());

    mccSpecies_.setSize(typeIds_.size());

    vibrationalETotal_.setSize(typeIds_.size());

    forAll(vibrationalETotal_, i)
    {
        vibrationalETotal_[i].setSize
        (
            cloud_.constProps(typeIds_[i]).vibrationalDoF()
        );
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

    nParcelsXnParticle_.setSize(typeIds_.size());

    for (auto& n : nParcelsXnParticle_)
    {
        n.setSize(mesh_.nCells());
    }

    mccSpecies_.setSize(typeIds_.size());

    for (auto& m : mccSpecies_)
    {
        m.setSize(mesh_.nCells());
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ShakhovRelaxation::calculateProperties()
{

    auto& cm = cloud_.cellPropMeasurements();

    forAll(typeIds_, iD)
    {

        rhoNMean_ = cm.rhoNMean()[iD];
        rhoMMean_ = cm.rhoMMean()[iD];
        linearKEMean_ = cm.linearKEMean()[iD];
        momentumMean_ = cm.momentumMean()[iD];
        rotationalEMean_ = cm.rotationalEMean()[iD];
        rotationalDofMean_ = cm.rotationalDofMean()[iD];
        electronicETotal_[iD] = cm.electronicETotal()[iD];
        rhoNMeanXnParticle_ = cm.rhoNMeanXnParticle()[iD];
        rhoMMeanXnParticle_ = cm.rhoMMeanXnParticle()[iD];
        momentumMeanXnParticle_ = cm.momentumMeanXnParticle()[iD];
        linearKEMeanXnParticle_ = cm.linearKEMeanXnParticle()[iD];

        muu_ = cm.muu()[iD];
        muv_ = cm.muv()[iD];
        muw_ = cm.muw()[iD];
        mvv_ = cm.mvv()[iD];
        mvw_ = cm.mvw()[iD];
        mww_ = cm.mww()[iD];
        mcc_ = cm.mcc()[iD];
        mccu_ = cm.mccu()[iD];
        mccv_ = cm.mccv()[iD];
        mccw_ = cm.mccw()[iD];

        eu_ = cm.eu()[iD];
        ev_ = cm.ev()[iD];
        ew_ = cm.ew()[iD];
        e_ = cm.e()[iD];

        rhoNMeanInt_ = cm.rhoNMeanInt()[iD];
        molsElec_ = cm.molsElec()[iD];

        nParcels_[iD] = cm.nParcels()[iD];
        nParcelsXnParticle_[iD] = cm.nParcelsXnParticle()[iD];
        mccSpecies_[iD] = cm.mccSpecies()[iD];

        nGroundElectronicLevel_[iD] = cm.nGroundElectronicLevel()[iD];
        nFirstElectronicLevel_[iD] = cm.nFirstElectronicLevel()[iD];

        forAll(vibrationalETotal_[iD], v)
        {
            vibrationalETotal_[iD][v] = cm.vibrationalETotal()[iD][v];
        }

    }

    forAll(mesh_.cells(), cell)
    {
        if (rhoNMean_[cell] > VSMALL)
        {

            const scalar cellVolume = mesh_.cellVolumes()[cell];

            // Number density
            rhoN_[cell] = rhoNMeanXnParticle_[cell]/cellVolume;

            // Velocity
            scalar rhoMMean = rhoMMeanXnParticle_[cell]/cellVolume;
            UMean_[cell] = momentumMeanXnParticle_[cell]/(rhoMMean*cellVolume);

            // Translational temperature
            scalar linearKEMean = 0.5*linearKEMeanXnParticle_[cell]/cellVolume;
            scalar rhoNMean = rhoNMeanXnParticle_[cell]/cellVolume;

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

        }
        else
        {
            rhoN_[cell] = 0.0;
            UMean_[cell] = vector::zero;
            p_[cell] = 0.0;
            translationalT_[cell] = 0.0;
            heatFluxVector_[cell] = vector::zero;
            pressureTensor_[cell] = tensor::zero;
        }

        // Rotational temperature
        if (rotationalDofMean_[cell] > VSMALL)
        {
            rotationalT_[cell] = (2.0/physicoChemical::k.value())*(rotationalEMean_[cell]/rotationalDofMean_[cell]);
        }
        else
        {
            rotationalT_[cell] = 0.0;
        }

         // Vibrational temperature

        scalarList degreesOfFreedomSpecies(typeIds_.size(), 0.0);
        scalarList vibTID(vibrationalETotal_.size(), 0.0);
        scalarField vibT(mesh_.nCells(), scalar(0.0));
        scalarField totalvDof(mesh_.nCells(), scalar(0.0));
        scalarField totalvDofOverall(mesh_.nCells(), scalar(0.0));

        List<scalarList> dofMode;
        List<scalarList> vibTMode;

        dofMode.setSize(typeIds_.size());
        vibTMode.setSize(typeIds_.size());

        forAll(dofMode, iD)
        {
            dofMode[iD].setSize
            (
                cloud_.constProps(typeIds_[iD]).vibrationalDoF(),
                0.0
            );

            vibTMode[iD].setSize
            (
                cloud_.constProps(typeIds_[iD]).vibrationalDoF(),
                0.0
            );
        }

        forAll(vibrationalETotal_, iD)
        {
            forAll(vibrationalETotal_[iD], v)
            {
                if
                (
                    vibrationalETotal_[iD][v][cell] > VSMALL
                 && nParcels_[iD][cell] > VSMALL
                 && dofMode.size() > VSMALL
                )
                {
                    scalar thetaV =
                        cloud_.constProps(typeIds_[iD]).thetaV()[v];

                    scalar vibrationalEMean =
                        vibrationalETotal_[iD][v][cell]
                       /nParcels_[iD][cell];

                    scalar iMean =
                        vibrationalEMean
                       /(physicoChemical::k.value()*thetaV);

                    vibTMode[iD][v] = thetaV / log(1.0 + (1.0/iMean));

                    dofMode[iD][v] =
                        (2.0*thetaV/vibTMode[iD][v])
                       /(exp(thetaV/vibTMode[iD][v]) - 1.0);
                }
            }

            forAll(dofMode[iD], v)
            {
                degreesOfFreedomSpecies[iD] += dofMode[iD][v];
            }

            forAll(dofMode[iD], v)
            {
                if (degreesOfFreedomSpecies[iD] > VSMALL)
                {
                    vibTID[iD] +=
                        vibTMode[iD][v]
                       *dofMode[iD][v]
                       /degreesOfFreedomSpecies[iD];
                }
            }


            totalvDof[cell] += degreesOfFreedomSpecies[iD];

            if
            (
                rhoNMeanInt_[cell] > VSMALL
             && rhoNMean_[cell] > VSMALL
             && nParcels_[iD][cell] > VSMALL
            )
            {
                scalar fraction =
                    nParcels_[iD][cell]
                   /rhoNMeanInt_[cell];

                scalar fractionOverall =
                    nParcels_[iD][cell]
                   /rhoNMean_[cell];

                totalvDofOverall[cell] +=
                    totalvDof[cell]
                   *(fractionOverall/fraction);

                vibT[cell] += vibTID[iD]*fraction;
            }
        }

        vibrationalT_[cell] = vibT[cell];

        // electronic temperature
        scalar totalEDof = 0.0;
        scalar elecT = 0.0;

        forAll(nParcels_, iD)
        {
            const scalarList& electronicEnergies =
                cloud_.constProps(typeIds_[iD]).electronicEnergyList();
            const labelList& degeneracies =
                cloud_.constProps(typeIds_[iD]).degeneracyList();

            if
            (
                nGroundElectronicLevel_[iD][cell] > VSMALL
             && nFirstElectronicLevel_[iD][cell] > VSMALL
             && nFirstElectronicLevel_[iD][cell]*degeneracies[0] !=
                nGroundElectronicLevel_[iD][cell]*degeneracies[1]
            )
            {

                scalar elecTID =
                    (electronicEnergies[1]-electronicEnergies[0])/
                    (
                        physicoChemical::k.value()*
                        log((nGroundElectronicLevel_[iD][cell]*
                         degeneracies[1])/
                        (nFirstElectronicLevel_[iD][cell]*
                        degeneracies[0]))
                    );


                scalar fraction = nParcels_[iD][cell]/molsElec_[cell];

                if (elecTID > VSMALL)
                {
                    elecT += fraction*elecTID;
                }


                scalar eDof =
                    (
                        2.0*(electronicETotal_[iD][cell]
                       /nParcels_[iD][cell])
                    )
                   /(physicoChemical::k.value()*elecTID);

                totalEDof += fraction*eDof;
            }
        }

        electronicT_[cell] = elecT;

        scalar nRotDof = 0.0;

        if (rhoNMean_[cell] > VSMALL)
        {
            nRotDof = rotationalDofMean_[cell] / rhoNMean_[cell];
        }

        overallT_[cell] =
            (
                (3.0*translationalT_[cell])
              + (nRotDof*rotationalT_[cell])
              + (totalvDof[cell]*vibrationalT_[cell])
              + (totalEDof*electronicT_[cell])
            )
           /(3.0 + nRotDof + totalvDof[cell] + totalEDof);

        // Relaxation frequency !!!Check mixtures and vibrational-electronic DoF
        Prandtl_[cell] = 0.0;
        scalar viscosity = 0.0;
        List<scalar> speciesVisc(typeIds_.size(), Zero);
        List<scalar> speciesPrandtl(typeIds_.size(), Zero);

        if (translationalT_[cell] > VSMALL)
        {
            forAll(typeIds_, iD)
            {

                const scalar& Tref = cloud_.constProps(iD).Tref();
                const scalar& mass = cloud_.constProps(iD).mass();
                const scalar& omega = cloud_.constProps(iD).omega();
                const scalar& d = cloud_.constProps(iD).d();
                const scalar& rotDoF = cloud_.constProps(iD).rotationalDoF();

                scalar speciesViscRef = 
                    7.5*sqrt(mass*physicoChemical::k.value()*Tref)
                    /(sqrt(mathematical::pi)*(5.0-2.0*omega)*(7.0-2.0*omega)*sqr(d));
                speciesVisc[iD] = speciesViscRef*pow(translationalT_[cell]/Tref,omega);
                viscosity += nParcels_[iD][cell]*speciesVisc[iD];

                speciesPrandtl[iD] += (5+rotDoF)/(7.5+rotDoF);
                Prandtl_[cell] += nParcels_[iD][cell]*speciesPrandtl[iD];

            }
            viscosity /= rhoNMean_[cell];
            Prandtl_[cell] /= rhoNMean_[cell]; 

            relaxFreq_[cell] = p_[cell]/viscosity;
        }
        else
        {
            Prandtl_[cell] = 0.0;
            relaxFreq_[cell] = 0.0;
        }

    }

    //Correct boundary conditions
    Prandtl_.correctBoundaryConditions();
    relaxFreq_.correctBoundaryConditions();
    rhoN_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    translationalT_.correctBoundaryConditions();
    rotationalT_.correctBoundaryConditions();
    vibrationalT_.correctBoundaryConditions();
    electronicT_.correctBoundaryConditions();
    overallT_.correctBoundaryConditions();
    UMean_.correctBoundaryConditions();
    pressureTensor_.correctBoundaryConditions();

}

void Foam::ShakhovRelaxation::resetProperties()
{

    forAll(mesh_.cells(), cell)
    {

        rhoMMean_[cell] = 0.0;
        linearKEMean_[cell] = 0.0;
        momentumMean_[cell] = vector::zero;
        rotationalEMean_[cell] = 0.0;
        rotationalDofMean_[cell] = 0.0;
        rhoNMeanXnParticle_[cell] = 0.0;
        rhoMMeanXnParticle_[cell] = 0.0;
        momentumMeanXnParticle_[cell] = vector::zero;
        linearKEMeanXnParticle_[cell] = 0.0;

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

        rhoNMeanInt_[cell] = 0.0;
        molsElec_[cell] = 0.0;

        forAll(typeIds_, iD)
        {

            nParcels_[iD][cell] = 0.0;
            nParcelsXnParticle_[iD][cell] = 0.0;
            mccSpecies_[iD][cell] = 0.0;
            electronicETotal_[iD][cell] = 0.0;
            nGroundElectronicLevel_[iD][cell] = 0.0;
            nFirstElectronicLevel_[iD][cell] = 0.0;

            forAll(vibrationalETotal_[iD], v)
            {
                vibrationalETotal_[iD][v][cell] = 0.0;
            }
        
        }

    }

}

void Foam::ShakhovRelaxation::relax()
{

    const scalar& deltaT = cloud_.mesh().time().deltaTValue();
    const scalar& nParticle = cloud_.nParticle();
    
    // Relax particles to local distribution
    label relaxations = 0;

    List<DynamicList<uspParcel*>>&
        cellOccupancy = cloud_.cellOccupancy();

    // Calculate required macroscopic properties
    calculateProperties();

    // Create macroscopic quantities interpolations
    //dictionary interpolationDict = mesh_.schemesDict().subDict("interpolationSchemes");;
    //autoPtr <Foam::interpolation<scalar>> PrandtlInterp = Foam::interpolationCellPoint<scalar>::New(interpolationDict, Prandtl_);
    //autoPtr <Foam::interpolation<scalar>> pInterp = Foam::interpolationCellPoint<scalar>::New(interpolationDict, p_);
    //autoPtr <Foam::interpolation<scalar>> translationalTInterp = Foam::interpolationCellPoint<scalar>::New(interpolationDict, translationalT_);
    //autoPtr <Foam::interpolation<vector>> UMeanInterp = Foam::interpolationCellPoint<vector>::New(interpolationDict, UMean_);
    //autoPtr <Foam::interpolation<vector>> heatFluxVectorInterp = Foam::interpolationCellPoint<vector>::New(interpolationDict, heatFluxVector_);
    //autoPtr <Foam::interpolation<tensor>> pressureTensorInterp = Foam::interpolationCellPoint<tensor>::New(interpolationDict, pressureTensor_);

    forAll(cellOccupancy, cell)
    {

        if (cloud_.cellCollModel(cell) == cloud_.relCollModel() && cellOccupancy[cell].size() >= 2)
        {

            const DynamicList<uspParcel*>& cellParcels(cellOccupancy[cell]);

            forAll(cellParcels, i)
            {

                if (cloud_.rndGen().sample01<scalar>() < 1.0-exp(-relaxFreq_[cell]*deltaT)) //cloud_.rndGen().sample01<scalar>()
                {

                    uspParcel& parcel = *cellParcels[i];
                    const scalar& mass = cloud_.constProps(parcel.typeId()).mass();
                    const vector& position = parcel.position();

                    // Interpolate macroscopic properties
                    //scalar Prandtl = PrandtlInterp().interpolate(position, cell);
                    //scalar p = pInterp().interpolate(position, cell);
                    //scalar translationalT = translationalTInterp().interpolate(position, cell);
                    //vector UMean = UMeanInterp().interpolate(position, cell);
                    //vector heatFluxVector = heatFluxVectorInterp().interpolate(position, cell);
                    //tensor pressureTensor = pressureTensorInterp().interpolate(position, cell);

                    scalar Prandtl = Prandtl_[cell];
                    scalar p = p_[cell];
                    scalar translationalT = translationalT_[cell];
                    vector UMean = UMean_[cell];
                    vector heatFluxVector = heatFluxVector_[cell];
                    tensor pressureTensor = pressureTensor_[cell];

                    // Relax particle
                    parcel.U() = samplePostRelaxationVelocity
                            (
                                mass,
                                Prandtl,
                                p,
                                translationalT,
                                UMean,
                                heatFluxVector
                            );

                    relaxations++;
                
                }
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

void Foam::ShakhovRelaxation::conserveMomentumAndEnergy
(
    const label& cell
)
{

    // Energy conservation scheme
    const scalar nParticle = cloud_.nParticle();
    //const scalar cellVolume = mesh_.cellVolumes()[cell];
    const DynamicList<uspParcel*>& cellParcels(cloud_.cellOccupancy()[cell]);

    scalar linearKEMeanXnParticle = 0.0;
    vector momentumMeanXnParticle = vector::zero;

    forAll(cellParcels, i)
    {

        const uspParcel& p = *cellParcels[i];

        const scalar mass = cloud_.constProps(p.typeId()).mass();
        const vector& U = p.U();
        const scalar& CWF = cloud_.cellWF(p.cell());
        const scalar& RWF = cloud_.axiRWF(p.position());

        linearKEMeanXnParticle += mass*(U & U)*CWF*RWF*nParticle;
        momentumMeanXnParticle += mass*U*CWF*RWF*nParticle;

    }

    // Velocity
    vector postUMean = momentumMeanXnParticle/rhoMMeanXnParticle_[cell];

    // Translational temperature
    scalar postTranslationalT =
        1.0/(3.0*physicoChemical::k.value()*rhoNMeanXnParticle_[cell])
       *(
            linearKEMeanXnParticle
          - rhoMMeanXnParticle_[cell]*(postUMean & postUMean)
        );

    forAll(cellParcels, i)
    {
            uspParcel& p = *cellParcels[i];
            p.U() = UMean_[cell] + (p.U()-postUMean)*sqrt(translationalT_[cell]/postTranslationalT);
    }

}

Foam::vector Foam::ShakhovRelaxation::samplePostRelaxationVelocity
(   
    const scalar& m,
    const scalar& Pr,
    const scalar& p,
    const scalar& T,
    const vector& U,
    const vector& q
)
{

    scalar u0(cloud_.maxwellianMostProbableSpeed(T,m));

    // Compute AJ parameter
    scalar B = 0.0;
    forAll(q, i) 
    {
        if (B < fabs(q[i])) 
        {
            B = fabs(q[i]);
        }
    }
    B = 2.0*(1.0-Pr)*B/(p*u0);

    vector v;
    scalar C;
    scalar A = 1.0 + 30.0*B;
    do {

        v = cloud_.rndGen().GaussNormal<vector>()/sqrt(2.0);
        C = 1.0+2.0*(1.0-Pr)/(p*u0)*(q.x()*v.x()+q.y()*v.y()+q.z()*v.z())*((sqr(v.x())+sqr(v.y())+sqr(v.z()))/2.5-1.0);

    } while(A*cloud_.rndGen().sample01<scalar>()>C);

    return U+u0*v;

}

// ************************************************************************* //

