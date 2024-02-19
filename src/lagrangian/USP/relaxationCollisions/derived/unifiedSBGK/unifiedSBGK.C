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

#include "unifiedSBGK.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(unifiedSBGK, 0);

addToRunTimeSelectionTable(relaxationCollisionModel, unifiedSBGK, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unifiedSBGK::unifiedSBGK
(
    const dictionary& dict,
    const polyMesh& mesh,
    uspCloud& cloud
)
:
    relaxationCollisionModel(dict, mesh, cloud),
    propertiesDict_(dict.subDict("collisionProperties")),
    Tref_(propertiesDict_.get<scalar>("Tref")),
    theta_(propertiesDict_.getOrDefault<scalar>("theta", 1.0)),
    macroInterpolation_(propertiesDict_.getOrDefault<bool>("macroInterpolation", false)),
    infoCounter_(0),
    shufflePasses_(5),
    maxProbResetValue_(0.999),
    performRelaxation_(mesh_.nCells(), true),
    maxProbReset_(mesh_.nCells(), true),
    maxProb_(mesh_.nCells(), 1.0),
    rhoNMean_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    linearKEMean_(mesh_.nCells(), 0.0),
    momentumMean_(mesh_.nCells(), vector::zero),
    rotationalEMean_(mesh_.nCells(), 0.0),
    rotationalDofMean_(mesh_.nCells(), 0.0),
    rhoNMeanXnParticle_(mesh_.nCells(), 0.0),
    rhoMMeanXnParticle_(mesh_.nCells(), 0.0),
    momentumMeanXnParticle_(mesh_.nCells(), vector::zero),
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
    rhoNMeanInt_(mesh_.nCells(), 0.0),
    molsElec_(mesh_.nCells(), 0.0),
    nParcels_(),
    nParcelsXnParticle_(),
    mccSpecies_(),
    nGroundElectronicLevel_(),
    nFirstElectronicLevel_(),
    electronicETotal_(),
    vibrationalETotal_(),
    Prandtl_
    (
        IOobject
        (
            "Prandtl",
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
            "relaxFreq",
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
    rotationalT_
    (
        IOobject
        (
            "rotationalT",
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
            "vibrationalT",
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
            "electronicT",
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
            "overallT",
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
        dimensionedVector(dimVelocity, vector::zero),
        zeroGradientFvPatchScalarField::typeName
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
    heatFluxVectorPrevious_
    (
        IOobject
        (
            "heatFluxVectorPrevious",
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
        dimensionedTensor(dimPressure, tensor::zero),
        zeroGradientFvPatchScalarField::typeName
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
        dimensionedTensor(dimPressure, tensor::zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    shearStressTensorPrevious_
    (
        IOobject
        (
            "shearStressTensorPrevious",
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

    electronicETotal_.setSize(typeIds_.size());

    for (auto& e : electronicETotal_)
    {
        e.setSize(mesh_.nCells());
    }

    vibrationalETotal_.setSize(typeIds_.size());

    forAll(vibrationalETotal_, i)
    {
        vibrationalETotal_[i].setSize
        (
            cloud_.constProps(typeIds_[i]).vibrationalDoF()
        );
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::unifiedSBGK::calculateProperties()
{

    auto& cm = cloud_.cellPropMeasurements();

    forAll(typeIds_, iD)
    {

        rhoNMean_ += cm.rhoNMean()[iD];
        rhoMMean_ += cm.rhoMMean()[iD];
        linearKEMean_ += cm.linearKEMean()[iD];
        momentumMean_ += cm.momentumMean()[iD];
        rotationalEMean_ += cm.rotationalEMean()[iD];
        rotationalDofMean_ += cm.rotationalDofMean()[iD];
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
        nParcelsXnParticle_[iD] += cm.nParcelsXnParticle()[iD];
        mccSpecies_[iD] += cm.mccSpecies()[iD];

        nGroundElectronicLevel_[iD] += cm.nGroundElectronicLevel()[iD];
        nFirstElectronicLevel_[iD] += cm.nFirstElectronicLevel()[iD];
        electronicETotal_[iD] += cm.electronicETotal()[iD];

        forAll(vibrationalETotal_[iD], v)
        {
            vibrationalETotal_[iD][v] += cm.vibrationalETotal()[iD][v];
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

            // Shear stress tensor
            scalar scalarPressure = (pressureTensor_[cell].xx() + pressureTensor_[cell].yy() + pressureTensor_[cell].zz())/3.0;

            shearStressTensor_[cell] = pressureTensor_[cell];
            shearStressTensor_[cell].xx() -= scalarPressure;
            shearStressTensor_[cell].yy() -= scalarPressure;
            shearStressTensor_[cell].zz() -= scalarPressure;

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

            // Scale macroscopic properties
            if (rhoNMean_[cell] > 2.0)
            {
                p_[cell] = rhoNMean_[cell]/(rhoNMean_[cell]-1.0)*p_[cell];
                translationalT_[cell] = rhoNMean_[cell]/(rhoNMean_[cell]-1.0)*translationalT_[cell];
                heatFluxVector_[cell] = sqr(rhoNMean_[cell])/(rhoNMean_[cell]-1.0)/(rhoNMean_[cell]-2.0)*heatFluxVector_[cell];
                shearStressTensor_[cell] = rhoNMean_[cell]/(rhoNMean_[cell]-1.0)*shearStressTensor_[cell];
            }
            else
            {
                performRelaxation_[cell] = false;
            }

        }
        else
        {
            performRelaxation_[cell] = false;
            rhoN_[cell] = 0.0;
            p_[cell] = 0.0;
            translationalT_[cell] = 0.0;
            UMean_[cell] = vector::zero;
            heatFluxVector_[cell] = vector::zero;
            pressureTensor_[cell] = tensor::zero;
            shearStressTensor_[cell] = tensor::zero;
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
                    1.25*(1.0+a)*(2.0+a)*sqrt(mass*physicoChemical::k.value()*Tref_)
                    /(a*(5.0-2.0*omega)*(7.0-2.0*omega)*sqrt(mathematical::pi)*sqr(d));
                    
                viscosity += nParcels_[iD][cell]*speciesViscRef*pow(translationalT_[cell]/Tref_,omega);

                Prandtl_[cell] += nParcels_[iD][cell]*(5.0+rotDoF)/(7.5+rotDoF);

            }
            viscosity /= rhoNMean_[cell];
            Prandtl_[cell] /= rhoNMean_[cell]; 

            relaxFreq_[cell] = p_[cell]/viscosity;
        }
        else
        {
            performRelaxation_[cell] = false;
            Prandtl_[cell] = 0.0;
            relaxFreq_[cell] = 0.0;
        }

    }

    // time-average and scale heat flux vector and pressure/shear-stress tensor
    const scalar& deltaT = cloud_.mesh().time().deltaTValue();

    forAll(mesh_.cells(), cell)
    {
        heatFluxVector_[cell] = theta_*heatFluxVector_[cell]/(1.0+0.5*Prandtl_[cell]*relaxFreq_[cell]*deltaT) + (1.0-theta_)*heatFluxVectorPrevious_[cell];
        shearStressTensor_[cell] = theta_*shearStressTensor_[cell]/(1.0+0.5*relaxFreq_[cell]*deltaT) + (1.0-theta_)*shearStressTensorPrevious_[cell];
        heatFluxVectorPrevious_[cell] = heatFluxVector_[cell];
        shearStressTensorPrevious_[cell] = shearStressTensor_[cell];
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
    heatFluxVector_.correctBoundaryConditions();
    pressureTensor_.correctBoundaryConditions();
    shearStressTensor_.correctBoundaryConditions();

}

void Foam::unifiedSBGK::resetProperties()
{

    forAll(mesh_.cells(), cell)
    {

        rhoNMean_[cell] = 0.0;
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
            nGroundElectronicLevel_[iD][cell] = 0.0;
            nFirstElectronicLevel_[iD][cell] = 0.0;
            electronicETotal_[iD][cell] = 0.0;

            forAll(vibrationalETotal_[iD], v)
            {
                vibrationalETotal_[iD][v][cell] = 0.0;
            }
        
        }

        if (maxProbReset_[cell])
        { 
            maxProb_[cell] *= maxProbResetValue_; 
        }
        maxProbReset_[cell] = true;

        performRelaxation_[cell] = true;

    }

}

void Foam::unifiedSBGK::relax()
{

    const scalar& deltaT = cloud_.mesh().time().deltaTValue();

    // Relax particles to local distribution
    label relaxations = 0;

    List<DynamicList<uspParcel*>>& cellOccupancy = cloud_.cellOccupancy();

    // Calculate required macroscopic properties
    calculateProperties();

    // Create macroscopic quantities interpolations
    autoPtr <Foam::interpolation<scalar>> PrandtlInterp;
    autoPtr <Foam::interpolation<scalar>> relaxFreqInterp;
    autoPtr <Foam::interpolation<scalar>> pInterp;
    autoPtr <Foam::interpolation<scalar>> translationalTInterp;
    autoPtr <Foam::interpolation<vector>> UMeanInterp;
    autoPtr <Foam::interpolation<vector>> heatFluxVectorInterp;
    autoPtr <Foam::interpolation<tensor>> shearStressTensorInterp;

    if (macroInterpolation_)
    {
        dictionary interpolationDict = mesh_.schemesDict().subDict("interpolationSchemes");
        PrandtlInterp = Foam::interpolationCellPoint<scalar>::New(interpolationDict, Prandtl_);
        relaxFreqInterp = Foam::interpolationCellPoint<scalar>::New(interpolationDict, relaxFreq_);
        pInterp = Foam::interpolationCellPoint<scalar>::New(interpolationDict, p_);
        translationalTInterp = Foam::interpolationCellPoint<scalar>::New(interpolationDict, translationalT_);
        UMeanInterp = Foam::interpolationCellPoint<vector>::New(interpolationDict, UMean_);
        heatFluxVectorInterp = Foam::interpolationCellPoint<vector>::New(interpolationDict, heatFluxVector_);
        shearStressTensorInterp = Foam::interpolationCellPoint<tensor>::New(interpolationDict, shearStressTensor_);
    }

    forAll(cellOccupancy, cell)
    {

        if (cloud_.cellCollModel(cell) == cloud_.relCollModel() && performRelaxation_[cell])
        {

            DynamicList<uspParcel*>& cellParcels(cellOccupancy[cell]);

            for (label pass = 1; pass <= shufflePasses_; pass++)
            {
                cloud_.rndGen().shuffle(cellParcels);
            }

            scalar parcelRelaxations = rhoNMean_[cell]*(1.0-exp(-relaxFreq_[cell]*deltaT));
            label nParcelRelaxations = int(parcelRelaxations);
            if (cloud_.rndGen().sample01<scalar>() < (parcelRelaxations-nParcelRelaxations))
            {
                nParcelRelaxations++;
            }
            nParcelRelaxations = min(nParcelRelaxations,rhoNMean_[cell]);

            for(label i = 0; i < nParcelRelaxations; i++)
            {

                uspParcel& parcel = *cellParcels[i];
                const scalar& mass = cloud_.constProps(parcel.typeId()).mass();
                const vector& position = parcel.position();

                // Relax particle
                if (macroInterpolation_)
                {
                    parcel.U() = samplePostRelaxationVelocity
                            (
                                cell,
                                mass,
                                PrandtlInterp().interpolate(position, cell),
                                relaxFreqInterp().interpolate(position, cell),
                                pInterp().interpolate(position, cell),
                                translationalTInterp().interpolate(position, cell),
                                UMeanInterp().interpolate(position, cell),
                                heatFluxVectorInterp().interpolate(position, cell),
                                shearStressTensorInterp().interpolate(position, cell)
                            );
                }
                else
                {
                    parcel.U() = samplePostRelaxationVelocity
                            (
                                cell,
                                mass,
                                Prandtl_[cell],
                                relaxFreq_[cell],
                                p_[cell],
                                translationalT_[cell],
                                UMean_[cell],
                                heatFluxVector_[cell],
                                shearStressTensor_[cell]
                            );                    
                }

                relaxations++;
                cloud_.cellPropMeasurements().nColls()[cell]++;
                
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

void Foam::unifiedSBGK::conserveMomentumAndEnergy
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
        rhoNMean_[cell]/(3.0*(rhoNMean_[cell]-1.0)*physicoChemical::k.value()*rhoNMeanXnParticle_[cell])
       *(
            linearKEMeanXnParticle
          - rhoMMeanXnParticle_[cell]*(postUMean & postUMean)
        );

    if (postTranslationalT > VSMALL)
    {
        forAll(cellParcels, i)
        {
                uspParcel& p = *cellParcels[i];
                p.U() = UMean_[cell] + (p.U()-postUMean)*sqrt(translationalT_[cell]/postTranslationalT);
        }
    }

}

Foam::vector Foam::unifiedSBGK::samplePostRelaxationVelocity
(   
    const label& cell,
    const scalar& m,
    const scalar& Pr,
    const scalar& rF,
    const scalar& p,
    const scalar& T,
    const vector& U,
    const vector& q,
    const tensor& s
)
{

    scalar u0(cloud_.maxwellianMostProbableSpeed(T,m));
    scalar tau = 0.5*rF*cloud_.mesh().time().deltaTValue();
    //scalar coeffQ = 2.0*(1.0-Pr*tau*(1.0+2.0/(exp(2.0*tau)-1.0)));
    //scalar coeffS = (1.0-tau*(1.0+2.0/(exp(2.0*tau)-1.0)));

    scalar factor = exp(-0.1/(2.0*tau));
    scalar coeffQ = 2.0*(factor*(1.0-Pr*tau*(1.0+2.0/(exp(2.0*tau)-1.0))) + (1.0-factor)*(exp(2.0*(1.0-Pr)*tau)-1.0)/(exp(2.0*tau)-1.0));
    scalar coeffS = factor*(1.0-tau*(1.0+2.0/(exp(2.0*tau)-1.0)));

    vector v;
    scalar prob;
    scalar vSq;
    scalar vTr;
    while(true)
    {

        v = cloud_.rndGen().GaussNormal<vector>()/sqrt(2.0);
        
        vSq = sqr(v.x())+sqr(v.y())+sqr(v.z());
        vTr = vSq/3.0;
        prob = 1.0+coeffQ/(p*u0)*(q.x()*v.x()+q.y()*v.y()+q.z()*v.z())*(vSq/2.5-1.0)
               +coeffS/p*(s.xx()*(v.x()*v.x()-vTr)+s.yy()*(v.y()*v.y()-vTr)+s.zz()*(v.z()*v.z()-vTr)+2.0*s.xy()*v.x()*v.y()+2.0*s.xz()*v.x()*v.z()+2.0*s.yz()*v.y()*v.z());

        if (prob > maxProb_[cell])
        {
            maxProb_[cell] = prob;
            maxProbReset_[cell] = false;
            break;
        }
        if (cloud_.rndGen().sample01<scalar>() < prob/maxProb_[cell])
        {
            break;
        }

    }

    return U+u0*v;

}

const Foam::dictionary&
Foam::unifiedSBGK::propertiesDict() const
{
    return propertiesDict_;
}

// ************************************************************************* //

