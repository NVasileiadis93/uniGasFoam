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

#include "EllipsoidalStatisticalRelaxation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(EllipsoidalStatisticalRelaxation, 0);

addToRunTimeSelectionTable(relaxationModel, EllipsoidalStatisticalRelaxation, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EllipsoidalStatisticalRelaxation::EllipsoidalStatisticalRelaxation
(
    const dictionary& dict,
    const polyMesh& mesh,
    uspCloud& cloud
)
:
    relaxationModel(dict, mesh, cloud),
    infoCounter_(0),
    rhoNMean_(),
    rhoNMeanXnParticle_(),
    rhoNMeanInt_(),
    molsElec_(),
    rhoMMean_(),
    rhoMMeanXnParticle_(),
    linearKEMean_(),
    linearKEMeanXnParticle_(),
    rotationalEMean_(),
    rotationalDofMean_(),
    muu_(),
    muv_(),
    muw_(),
    mvv_(),
    mvw_(),
    mww_(),
    mcc_(),
    mccu_(),
    mccv_(),
    mccw_(),
    eu_(),
    ev_(),
    ew_(),
    e_(),
    momentumMean_(),
    momentumMeanXnParticle_(),
    vibrationalETotal_(),
    electronicETotal_(),
    nParcels_(),
    nParcelsXnParticle_(),
    mccSpecies_(),
    nGroundElectronicLevel_(),
    nFirstElectronicLevel_(),
    viscosity_(),
    Prandtl_(),
    relaxFreq_(),
    rhoN_(),
    p_(),
    translationalT_(),
    rotationalT_(),
    vibrationalT_(),
    electronicT_(),
    overallT_(),
    UMean_(),
    pressureTensor_()
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

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::EllipsoidalStatisticalRelaxation::calculateProperties
(
    const label& cell
)
{

    auto& cm = cloud_.cellPropMeasurements();

    forAll(typeIds_, iD)
    {

        rhoNMean_ = cm.rhoNMean()[iD][cell];
        rhoMMean_ = cm.rhoMMean()[iD][cell];
        linearKEMean_ = cm.linearKEMean()[iD][cell];
        momentumMean_ = cm.momentumMean()[iD][cell];
        rotationalEMean_ = cm.rotationalEMean()[iD][cell];
        rotationalDofMean_ = cm.rotationalDofMean()[iD][cell];
        electronicETotal_[iD] = cm.electronicETotal()[iD][cell];
        rhoNMeanXnParticle_ = cm.rhoNMeanXnParticle()[iD][cell];
        rhoMMeanXnParticle_ = cm.rhoMMeanXnParticle()[iD][cell];
        momentumMeanXnParticle_ = cm.momentumMeanXnParticle()[iD][cell];
        linearKEMeanXnParticle_ = cm.linearKEMeanXnParticle()[iD][cell];

        muu_ = cm.muu()[iD][cell];
        muv_ = cm.muv()[iD][cell];
        muw_ = cm.muw()[iD][cell];
        mvv_ = cm.mvv()[iD][cell];
        mvw_ = cm.mvw()[iD][cell];
        mww_ = cm.mww()[iD][cell];
        mcc_ = cm.mcc()[iD][cell];
        mccu_ = cm.mccu()[iD][cell];
        mccv_ = cm.mccv()[iD][cell];
        mccw_ = cm.mccw()[iD][cell];

        eu_ = cm.eu()[iD][cell];
        ev_ = cm.ev()[iD][cell];
        ew_ = cm.ew()[iD][cell];
        e_ = cm.e()[iD][cell];

        rhoNMeanInt_ = cm.rhoNMeanInt()[iD][cell];
        molsElec_ = cm.molsElec()[iD][cell];

        nParcels_[iD] = cm.nParcels()[iD][cell];
        nParcelsXnParticle_[iD] = cm.nParcelsXnParticle()[iD][cell];
        mccSpecies_[iD] = cm.mccSpecies()[iD][cell];

        nGroundElectronicLevel_[iD] = cm.nGroundElectronicLevel()[iD][cell];
        nFirstElectronicLevel_[iD] = cm.nFirstElectronicLevel()[iD][cell];

        forAll(vibrationalETotal_[iD], v)
        {
            vibrationalETotal_[iD][v] = cm.vibrationalETotal()[iD][v][cell];
        }

    }

    if (rhoNMean_ > VSMALL)
    {
        const scalar cellVolume = mesh_.cellVolumes()[cell];

        // Number density
        rhoN_ = rhoNMeanXnParticle_/cellVolume;


        // Velocity
        scalar rhoMMean = rhoMMeanXnParticle_/cellVolume;
        UMean_ = momentumMeanXnParticle_/(rhoMMean*cellVolume);

        // Translational temperature
        scalar linearKEMean = 0.5*linearKEMeanXnParticle_/cellVolume;
        scalar rhoNMean = rhoNMeanXnParticle_/cellVolume;

        translationalT_ =
            2.0/(3.0*physicoChemical::k.value()*rhoNMean)
           *(
                linearKEMean
              - 0.5*rhoMMean*(UMean_ & UMean_)
            );

        // Pressure
        p_ =
            rhoN_*physicoChemical::k.value()
            *translationalT_;

        // Pressure tensor
        pressureTensor_.xx() = rhoN_*
        (
            muu_/(rhoNMean_) -
            (
                (rhoMMean_/(rhoNMean_))
                *UMean_.x()*UMean_.x()
            )
        );
        pressureTensor_.xy() = rhoN_*
        (
            muv_/(rhoNMean_) -
            ((rhoMMean_/(rhoNMean_)))
            *UMean_.x()*UMean_.y()

        );
        pressureTensor_.xz() = rhoN_*
        (
            muw_/(rhoNMean_) -
            ((rhoMMean_/(rhoNMean_))
            *UMean_.x()*UMean_.z())
        );
        pressureTensor_.yx() = pressureTensor_.xy();
        pressureTensor_.yy() = rhoN_*
        (
            mvv_/(rhoNMean_) -
            ((rhoMMean_/(rhoNMean_)))
            *UMean_.y()*UMean_.y()
        );
        pressureTensor_.yz() = rhoN_*
        (
            mvw_/(rhoNMean_) -
            ((rhoMMean_/(rhoNMean_))
            *UMean_.y()*UMean_.z())
        );
        pressureTensor_.zx() = pressureTensor_.xz();
        pressureTensor_.zy() = pressureTensor_.yz();
        pressureTensor_.zz() = rhoN_*
        (
            mww_/(rhoNMean_) -
            ((rhoMMean_/(rhoNMean_))
            *UMean_.z()*UMean_.z())
        );

    }
    else
    {
        rhoN_ = 0.0;
        UMean_ = vector::zero;
        p_ = 0.0;
        translationalT_ = 0.0;
        pressureTensor_ = tensor::zero;
    }

    // Rotational temperature
    if (rotationalDofMean_ > VSMALL)
    {
        rotationalT_ = (2.0/physicoChemical::k.value())*(rotationalEMean_/rotationalDofMean_);
    }
    else
    {
        rotationalT_ = 0.0;
    }

    // Vibrational temperature

    scalar vibT = 0.0;
    scalar totalvDof = 0.0;
    scalar totalvDofOverall = 0.0;
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
            if (vibrationalETotal_[iD][v] > VSMALL && nParcels_[iD] > VSMALL && dofMode.size() > VSMALL
            )
            {
                scalar thetaV = cloud_.constProps(typeIds_[iD]).thetaV()[v];

                scalar vibrationalEMean = vibrationalETotal_[iD][v]/nParcels_[iD];

                scalar iMean = vibrationalEMean/(physicoChemical::k.value()*thetaV);

                vibTMode[iD][v] = thetaV / log(1.0 + (1.0/iMean));

                dofMode[iD][v] = (2.0*thetaV/vibTMode[iD][v])/(exp(thetaV/vibTMode[iD][v]) - 1.0);
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
                vibTID[iD] += vibTMode[iD][v]*dofMode[iD][v]/degreesOfFreedomSpecies[iD];
            }
        }

        totalvDof += degreesOfFreedomSpecies[iD];

        if(rhoNMeanInt_ > VSMALL && rhoNMean_ > VSMALL && nParcels_[iD] > VSMALL)
        {
            scalar fraction = nParcels_[iD]/rhoNMeanInt_;

            scalar fractionOverall = nParcels_[iD]/rhoNMean_;

            totalvDofOverall += totalvDof*(fractionOverall/fraction);

            vibT += vibTID[iD]*fraction;
        }
    }

    vibrationalT_ = vibT;

    // Electronic temperature
    scalar totalEDof = 0.0;
    scalar elecT = 0.0;

    forAll(nParcels_, iD)
    {
        const scalarList& electronicEnergies = cloud_.constProps(typeIds_[iD]).electronicEnergyList();
        const labelList& degeneracies = cloud_.constProps(typeIds_[iD]).degeneracyList();

        if(
            nGroundElectronicLevel_[iD] > VSMALL
         && nFirstElectronicLevel_[iD] > VSMALL
         && nFirstElectronicLevel_[iD]*degeneracies[0] != nGroundElectronicLevel_[iD]*degeneracies[1]
        )
        {

            scalar elecTID =
                (electronicEnergies[1]-electronicEnergies[0])/
                (
                    physicoChemical::k.value()*
                    log((nGroundElectronicLevel_[iD]*
                     degeneracies[1])/
                    (nFirstElectronicLevel_[iD]*
                    degeneracies[0]))
                );


            scalar fraction = nParcels_[iD]/molsElec_;

            if (elecTID > VSMALL)
            {
                elecT += fraction*elecTID;
            }

            scalar eDof =
                (
                    2.0*(electronicETotal_[iD]
                   /nParcels_[iD])
                )
               /(physicoChemical::k.value()*elecTID);

            totalEDof += fraction*eDof;
        }
    }

    electronicT_ = elecT;

    scalar nRotDof = 0.0;

    if (rhoNMean_ > VSMALL)
    {
        nRotDof = rotationalDofMean_ / rhoNMean_;
    }

    // Overall temperature

    overallT_ =
        (
            (3.0*translationalT_)
          + (nRotDof*rotationalT_)
          + (totalvDof*vibrationalT_)
          + (totalEDof*electronicT_)
        )
       /(3.0 + nRotDof + totalvDof + totalEDof);

    

    // Relaxation frequency !!!Check mixtures and vibrational-electronic DoF

    viscosity_ = 0.0;
    Prandtl_ = 0.0;
    List<scalar> speciesVisc(typeIds_.size(), Zero);
    List<scalar> speciesPrandtl(typeIds_.size(), Zero);

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
        speciesVisc[iD] = speciesViscRef*pow(translationalT_/Tref,omega);
        viscosity_ += nParcels_[iD]*speciesVisc[iD];

        speciesPrandtl[iD] += (5+rotDoF)/(7.5+rotDoF);
        Prandtl_ += nParcels_[iD]*speciesPrandtl[iD];

    }
    
    viscosity_ /= rhoNMean_;
    Prandtl_ /= rhoNMean_; 

    //Prandtl_[cell] = 1.0; // reduce to BGK
    relaxFreq_ = Prandtl_*p_/viscosity_;

}

void Foam::EllipsoidalStatisticalRelaxation::resetProperties
(
    const label& cell
)
{

    rhoNMean_ = 0.0;
    rhoMMean_ = 0.0;
    linearKEMean_ = 0.0;
    momentumMean_ = vector::zero;
    rotationalEMean_ = 0.0;
    rotationalDofMean_ = 0.0;
    rhoNMeanXnParticle_ = 0.0;
    rhoMMeanXnParticle_ = 0.0;
    momentumMeanXnParticle_ = vector::zero;
    linearKEMeanXnParticle_ = 0.0;

    muu_ = 0.0;
    muv_ = 0.0;
    muw_ = 0.0;
    mvv_ = 0.0;
    mvw_ = 0.0;
    mww_ = 0.0;
    mcc_ = 0.0;
    mccu_ = 0.0;
    mccv_ = 0.0;
    mccw_ = 0.0;

    eu_ = 0.0;
    ev_ = 0.0;
    ew_ = 0.0;
    e_ = 0.0;

    rhoNMeanInt_ = 0.0;
    molsElec_ = 0.0;

    forAll(typeIds_, iD)
    {

        nParcels_[iD] = 0.0;
        nParcelsXnParticle_[iD] = 0.0;
        mccSpecies_[iD] = 0.0;
        electronicETotal_[iD] = 0.0;
        nGroundElectronicLevel_[iD] = 0.0;
        nFirstElectronicLevel_[iD] = 0.0;

        forAll(vibrationalETotal_[iD], v)
        {
            vibrationalETotal_[iD][v] = 0.0;
        }
    
    }
}

void Foam::EllipsoidalStatisticalRelaxation::relax()
{

    const scalar& deltaT = cloud_.mesh().time().deltaTValue();
    const scalar& nParticle = cloud_.nParticle();
    
    // Relax particles to local distribution
    label relaxations = 0;

    List<DynamicList<uspParcel*>>&
        cellOccupancy = cloud_.cellOccupancy();

    forAll(cellOccupancy, cell)
    {

        if (cloud_.cellCollModel(cell) == cloud_.relCollModel())
        {
            // Calculate required macroscopic properties
            calculateProperties(cell);

            const DynamicList<uspParcel*>& cellParcels(cellOccupancy[cell]);

            forAll(cellParcels, i)
            {

                if (cloud_.rndGen().sample01<scalar>() < 1.0-exp(-relaxFreq_*deltaT)) //cloud_.rndGen().sample01<scalar>()
                {

                    uspParcel& p = *cellParcels[i];
                    const scalar mass = cloud_.constProps(p.typeId()).mass();

                    // Relax particle
                    p.U() = samplePostRelaxationVelocity
                            (
                                mass,
                                Prandtl_,
                                p_,
                                translationalT_,
                                UMean_,
                                pressureTensor_
                            );



                    relaxations++;
                
                }
            }

            conserveMomentumAndEnergy(cell);

            resetProperties(cell);

        }
    }

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

void Foam::EllipsoidalStatisticalRelaxation::conserveMomentumAndEnergy
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

    scalar postTranslationalT;
    vector postUMean;

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
    postUMean = momentumMeanXnParticle/rhoMMeanXnParticle_;

    // Translational temperature
    postTranslationalT =
        1.0/(3.0*physicoChemical::k.value()*rhoNMeanXnParticle_)
       *(
            linearKEMeanXnParticle
          - rhoMMeanXnParticle_*(postUMean & postUMean)
        );

    forAll(cellParcels, i)
    {
            uspParcel& p = *cellParcels[i];
            p.U() = UMean_ + (p.U()-postUMean)*sqrt(translationalT_/postTranslationalT);
    }

    //std::cout << "-------------------------------- ESBGK --------------------------------" << std::endl;
    //std::cout << Prandtl_ << " " << viscosity_ << " " << 1.0-exp(-relaxFreq_*cloud_.mesh().time().deltaTValue()) << std::endl;
    //std::cout << translationalT_ << " " << postTranslationalT << std::endl;
    //std::cout << UMean_.x() << " " << UMean_.y() << " " << UMean_.z() << std::endl;
    //std::cout << postUMean.x() << " " << postUMean.y() << " " << postUMean.z() << std::endl;
    //std::cin.get();
}

Foam::vector Foam::EllipsoidalStatisticalRelaxation::samplePostRelaxationVelocity
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

