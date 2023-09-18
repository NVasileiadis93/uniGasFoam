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
    timeSteps_(0),
    nAvTimeSteps_(0),
    decomposeInterval_(dict.subDict("hybridProperties").get<label>("decomposeInterval")),
    bMax_(dict.subDict("hybridProperties").get<scalar>("bMax")),
    Tref_(dict.subDict("collisionCoeffs").get<scalar>("Tref")),
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
    nParcelsXnParticle_(),
    mccSpecies_(),
    nGroundElectronicLevel_(),
    nFirstElectronicLevel_(),
    breakdownParameter_
    (
        IOobject
        (
            "breakdownParameter",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
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
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimMass*pow(dimTime,-3), Zero)
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
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero)
    )
{

    typeIds_.setSize(cloud_.typeIdList().size());
    forAll(typeIds_,iD)
    {
        typeIds_[iD] = iD;
    }

    nGroundElectronicLevel_.setSize(typeIds_.size());
    for (auto& f : nGroundElectronicLevel_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    nFirstElectronicLevel_.setSize(typeIds_.size());
    for (auto& f : nFirstElectronicLevel_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    vibrationalETotal_.setSize(typeIds_.size());

    electronicETotal_.setSize(typeIds_.size());
    for (auto& f : electronicETotal_)
    {
        f.setSize(mesh_.nCells(), 0.0);
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

    mccSpecies_.setSize(typeIds_.size());
    for (auto& f : mccSpecies_)
    {
        f.setSize(mesh_.nCells());
    }

}


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uspHybridDecomposer::update()
{

    timeSteps_++;

    nAvTimeSteps_++;

    // get cell measurements
    auto& cm = cloud_.cellPropMeasurements();

    forAll(cm.rhoNMean(), i)
    {

        rhoNMean_+= cm.rhoNMean()[i];
        rhoMMean_ += cm.rhoMMean()[i];
        linearKEMean_ += cm.linearKEMean()[i];
        momentumMean_ += cm.momentumMean()[i];
        rotationalEMean_ += cm.rotationalEMean()[i];
        rotationalDofMean_ += cm.rotationalDofMean()[i];
        electronicETotal_[i] += cm.electronicETotal()[i];
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

        rhoNMeanInt_ += cm.rhoNMeanInt()[i];
        molsElec_ += cm.molsElec()[i];

        nParcels_[i] += cm.nParcels()[i];
        nParcelsXnParticle_[i] += cm.nParcelsXnParticle()[i];
        mccSpecies_[i] += cm.mccSpecies()[i];

        nGroundElectronicLevel_[i] += cm.nGroundElectronicLevel()[i];
        nFirstElectronicLevel_[i] += cm.nFirstElectronicLevel()[i];

        forAll(vibrationalETotal_[i], v)
        {
            vibrationalETotal_[i][v] += cm.vibrationalETotal()[i][v];
        }

    }

    if (timeSteps_ == decomposeInterval_)
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

                // Rotational temperature
                if (rotationalDofMean_[cell] > VSMALL)
                {
                    scalar rotationalEMean =
                        rotationalEMean_[cell]/nAvTimeSteps_;
                    scalar rotationalDofMean =
                        rotationalDofMean_[cell]/nAvTimeSteps_;

                    rotationalT_[cell] =
                        (2.0/physicoChemical::k.value())
                       *(rotationalEMean/rotationalDofMean);
                }
                else
                {
                    rotationalT_[cell] = 0.0;
                }

                // Vibrational temperature
                scalarField vibT(mesh_.nCells(), scalar(0.0));
                scalarField totalvDof(mesh_.nCells(), scalar(0.0));
                scalarField totalvDofOverall(mesh_.nCells(), scalar(0.0));
                scalarList degreesOfFreedomSpecies(typeIds_.size(), 0.0);
                scalarList vibTID(vibrationalETotal_.size(), 0.0);

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
            }
            else
            {
                rhoN_[cell] = 0.0;
                rhoM_[cell] = 0.0;
                p_[cell] = 0.0;
                translationalT_[cell] = 0.0;
                rotationalT_[cell] = 0.0;
                vibrationalT_[cell] = 0.0;
                electronicT_[cell] = 0.0;
                UMean_[cell] = vector::zero;
                heatFluxVector_[cell] = vector::zero;
                pressureTensor_[cell] = tensor::zero;
                shearStressTensor_[cell] = tensor::zero;
            }
        }


        // Calculate breakdown parameter
        forAll(mesh_.cells(), cell)
        {

            scalar massMean = 0.0;
            forAll(typeIds_,iD)
            {
                massMean += cloud_.constProps(iD).mass()*nParcels_[iD][cell];
            }
            massMean /= rhoNMean_[cell];

            scalar u0(
                cloud_.maxwellianMostProbableSpeed
                (
                    translationalT_[cell],
                    massMean
                )
            );

            scalar maxHeatFlux = VSMALL;
            forAll(heatFluxVector_[cell], i) 
            {
                if (maxHeatFlux < fabs(heatFluxVector_[cell][i])) 
                {
                    maxHeatFlux = fabs(heatFluxVector_[cell][i]);
                }
            }
            maxHeatFlux = 2.0*maxHeatFlux/(p_[cell]*u0);

            scalar maxShearStress = VSMALL;
            forAll(shearStressTensor_[cell],i) 
            {
                if (maxShearStress < fabs(shearStressTensor_[cell][i])) 
                {
                    maxShearStress = fabs(shearStressTensor_[cell][i]);
                }
            }
            maxShearStress = maxShearStress/p_[cell];

            breakdownParameter_[cell] = max(maxHeatFlux,maxShearStress);

        }

        // update domain decomposition
        forAll(mesh_.cells(), cell)
        {
            if (breakdownParameter_[cell] > bMax_)
            {
                cloud_.cellCollModel()[cell] = cloud_.binCollModel();
            }
            else
            {
                cloud_.cellCollModel()[cell] = cloud_.relCollModel();
            }
        }

        // reset
        forAll(rhoN_, cell)
        {

        }
        timeSteps_ = 0;

    }

}

}  // End namespace Foam

// ************************************************************************* //

