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

#include "ChapmanEnskogRelaxation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(ChapmanEnskogRelaxation, 0);

addToRunTimeSelectionTable(relaxationModel, ChapmanEnskogRelaxation, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ChapmanEnskogRelaxation::ChapmanEnskogRelaxation
(
    const dictionary& dict,
    const polyMesh& mesh,
    uspCloud& cloud
)
:
    relaxationModel(dict, mesh, cloud),
    timestepCounter_(0),
    heatFluxX_(0),
    heatFluxY_(0),
    shearStressXY_(0),
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
    heatFluxVector_(),
    pressureTensor_(),
    stressTensor_()
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

void Foam::ChapmanEnskogRelaxation::calculateProperties
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

    if (sum(nParcels_) > VSMALL)
    {

        const scalar cellVolume = mesh_.cellVolumes()[cell];

        rhoN_ = rhoNMeanXnParticle_/cellVolume;

        scalar rhoMMean = rhoMMeanXnParticle_/cellVolume;
        UMean_ = momentumMeanXnParticle_/(rhoMMean*cellVolume);

        scalar linearKEMean = 0.5*linearKEMeanXnParticle_/cellVolume;
        scalar rhoNMean = rhoNMeanXnParticle_/cellVolume;
        translationalT_ =
            2.0/(3.0*physicoChemical::k.value()*rhoNMean)
           *(
                linearKEMean
              - 0.5*rhoMMean*(UMean_ & UMean_)
            );

        p_ = rhoN_*physicoChemical::k.value()*translationalT_;

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

        scalar scalarPressure = (1.0/3.0)*
                                (pressureTensor_.xx() +
                                pressureTensor_.yy() +
                                pressureTensor_.zz());

        stressTensor_ = -pressureTensor_;
        stressTensor_.xx() += scalarPressure;
        stressTensor_.yy() += scalarPressure;
        stressTensor_.zz() += scalarPressure;

        heatFluxVector_.x() = rhoN_*
        (
            0.5*(mccu_/(rhoNMean_)) -
            0.5*(mcc_/(rhoNMean_))*
            UMean_.x() + eu_/(rhoNMean_) -
            (e_/(rhoNMean_))*UMean_.x()
        ) -
            pressureTensor_.xx()*UMean_.x() -
            pressureTensor_.xy()*UMean_.y() -
            pressureTensor_.xz()*UMean_.z();

        heatFluxVector_.y() = rhoN_*
        (
            0.5*(mccv_/(rhoNMean_)) -
            0.5*(mcc_/(rhoNMean_))*
            UMean_.y() + ev_/(rhoNMean_)-
            (e_/(rhoNMean_))*UMean_.y()
        ) -
            pressureTensor_.yx()*UMean_.x() -
            pressureTensor_.yy()*UMean_.y() -
            pressureTensor_.yz()*UMean_.z();

        heatFluxVector_.z() = rhoN_*
        (
            0.5*(mccw_/(rhoNMean_)) -
            0.5*(mcc_/(rhoNMean_))*
            UMean_.z() + ew_/(rhoNMean_) -
            (e_/(rhoNMean_))*UMean_.z()
        ) -
            pressureTensor_.zx()*UMean_.x() -
            pressureTensor_.zy()*UMean_.y() -
            pressureTensor_.zz()*UMean_.z();

    }
    else
    {
        rhoN_ = 0.0;
        p_ = 0.0;
        translationalT_ = 0.0;
        UMean_ = vector::zero;
        heatFluxVector_ = vector::zero;
        stressTensor_ = tensor::zero;
    }

    // Relaxation frequency !!!Check mixtures and vibrational-electronic DoF
    viscosity_ = 0.0;
    Prandtl_ = 0.0;
    List<scalar> speciesVisc(typeIds_.size(), Zero);
    List<scalar> speciesPrandtl(typeIds_.size(), Zero);

    if (translationalT_ > VSMALL)
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
            speciesVisc[iD] = speciesViscRef*pow(translationalT_/Tref,omega);
            viscosity_ += nParcels_[iD]*speciesVisc[iD];

            speciesPrandtl[iD] += (5+rotDoF)/(7.5+rotDoF);
            Prandtl_ += nParcels_[iD]*speciesPrandtl[iD];

        }
        viscosity_ /= rhoNMean_;
        Prandtl_ /= rhoNMean_; 

        relaxFreq_ = Prandtl_*p_/viscosity_;
    }
    else
    {
        relaxFreq_ = 0.0;
    }

}

void Foam::ChapmanEnskogRelaxation::resetProperties()
{

    rhoNMean_ = scalar(0.0);
    rhoMMean_ = scalar(0.0);
    linearKEMean_ = scalar(0.0);
    momentumMean_ = vector::zero;
    rotationalEMean_ = scalar(0.0);
    rotationalDofMean_ = scalar(0.0);
    rhoNMeanInt_ = scalar(0.0);
    molsElec_ = scalar(0.0),
    muu_ = scalar(0.0);
    muv_ = scalar(0.0);
    muw_ = scalar(0.0);
    mvv_ = scalar(0.0);
    mvw_ = scalar(0.0);
    mww_ = scalar(0.0);
    mcc_ = scalar(0.0);
    mccu_ = scalar(0.0);
    mccv_ = scalar(0.0);
    mccw_ = scalar(0.0);
    eu_ = scalar(0.0);
    ev_ = scalar(0.0);
    ew_ = scalar(0.0);
    e_ = scalar(0.0);
    rhoNMeanXnParticle_ = scalar(0.0);
    rhoMMeanXnParticle_ = scalar(0.0);
    momentumMeanXnParticle_ = vector::zero;
    linearKEMeanXnParticle_ = scalar(0.0);

    forAll(electronicETotal_, iD)
    {

        electronicETotal_[iD] = 0.0;
        mccSpecies_[iD] = 0.0;
        nParcels_[iD] = 0.0;
        nGroundElectronicLevel_[iD] = 0.0;
        nFirstElectronicLevel_[iD] = 0.0;
        nParcelsXnParticle_[iD] = 0.0;

        forAll(vibrationalETotal_[iD], v)
        {
           vibrationalETotal_[iD][v] = 0.0;
        }

    }

}

void Foam::ChapmanEnskogRelaxation::relax()
{

    timestepCounter_++;

    const scalar& deltaT = cloud_.mesh().time().deltaTValue();
    const scalar& nParticle = cloud_.nParticle();
    
    // Relax particles to local distribution
    label relaxations = 0;

    List<DynamicList<uspParcel*>>&
        cellOccupancy = cloud_.cellOccupancy();

    std::ifstream macrofile;
    macrofile.open("./Couette_Quantities.dat");

    forAll(cellOccupancy, cell)
    {

        if (cloud_.cellCollModel(cell) == cloud_.relCollModel() && cellOccupancy[cell].size() >= 2)
        {

            // Calculate required macroscopic properties
            calculateProperties(cell);

            if (cell == 0)
            {
                heatFluxX_ += heatFluxVector_.x();
                heatFluxY_ += heatFluxVector_.y();
                shearStressXY_ += stressTensor_.xy();

                //pressureTensor_.xx() = rhoN_*
                //(
                //    muu_/(rhoNMean_) -
                //    (
                //        (rhoMMean_/(rhoNMean_))
                //        *UMean_.x()*UMean_.x()
                //    )
                //);

                std::cout << muv_ << " " << rhoNMean_ << " " << rhoMMean_ << " " << UMean_.x() << " " << UMean_.y() << " " << pressureTensor_.xy() << std::endl;

            }

            heatFluxVector_ = vector::zero;
            stressTensor_ = tensor::zero;

            //translationalT_ = 200;
            //UMean_ = vector(50, 100, -50);
            //heatFluxVector_ = vector(50, -20, -30);
            //stressTensor_ = tensor(1, 2, 3,
            //                       2, 2, 1,
            //                       3, 1, -3);
        
            //macrofile >> heatFluxVector_.x() >> heatFluxVector_.y() >> heatFluxVector_.z()
            //          >> stressTensor_.xx() >> stressTensor_.xy() >> stressTensor_.xz()
            //          >> stressTensor_.yx() >> stressTensor_.yy() >> stressTensor_.yz()
            //          >> stressTensor_.zx() >> stressTensor_.zy() >> stressTensor_.zz();

            const DynamicList<uspParcel*>& cellParcels(cellOccupancy[cell]);

            forAll(cellParcels, i)
            {

                //if (cloud_.rndGen().sample01<scalar>() < relaxFreq_*deltaT) //cloud_.rndGen().sample01<scalar>()
                //{

                    uspParcel& p = *cellParcels[i];
                    const scalar mass = cloud_.constProps(p.typeId()).mass();

                    scalar u0(
                        cloud_.maxwellianMostProbableSpeed
                        (
                            translationalT_,
                            mass
                        )
                    );

                    scalar maxHeatFlux = VSMALL;
                    forAll(heatFluxVector_,i) 
                    {
                        if (maxHeatFlux < fabs(heatFluxVector_[i])) 
                        {
                            maxHeatFlux = fabs(heatFluxVector_[i]);
                        }
                    }
                    maxHeatFlux = 2.0*maxHeatFlux/(p_*u0);

                    scalar maxStress = VSMALL;
                    forAll(stressTensor_,i) 
                    {
                        if (maxStress < fabs(stressTensor_[i])) 
                        {
                            maxStress = fabs(stressTensor_[i]);
                        }
                    }
                    maxStress = maxStress/p_;

                    scalar breakdownParameter = max(maxHeatFlux,maxStress);

                    // Relax particle
                    p.U() = samplePostRelaxationVelocity
                            (
                                breakdownParameter,
                                mass,
                                p_,
                                translationalT_,
                                UMean_,
                                heatFluxVector_,
                                stressTensor_
                            );

                    /*p.U() = samplePostRelaxationVelocityESBGK
                            (
                                mass,
                                Prandtl_,
                                p_,
                                translationalT_,
                                UMean_,
                                pressureTensor_
                            );*/

                    relaxations++;

                //}

            }
            
            conserveMomentumAndEnergy(cell);

            resetProperties();

        }
    }

    macrofile.close();

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

        std::ofstream resultfile;
        resultfile.open("./results.dat",std::ios_base::app);
        resultfile << heatFluxX_/scalar(timestepCounter_) << " "  << heatFluxY_/scalar(timestepCounter_) << " " << shearStressXY_/scalar(timestepCounter_) << std::endl;
        resultfile.close();
    }

    if (timestepCounter_ > 100000)
    {
        timestepCounter_ = 0;
        heatFluxX_ = 0.0;
        heatFluxY_ = 0.0;
        shearStressXY_ = 0.0;
    }

}

void Foam::ChapmanEnskogRelaxation::conserveMomentumAndEnergy
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

}

Foam::vector Foam::ChapmanEnskogRelaxation::samplePostRelaxationVelocity
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
    scalar A = 1.0 + 60.0*B;
    scalar u0(cloud_.maxwellianMostProbableSpeed(T,m));

    do {

        v = cloud_.rndGen().GaussNormal<vector>()/sqrt(2.0);
        gamma = 1.0+2.0/(p*u0)*(q.x()*v.x()+q.y()*v.y()+q.z()*v.z())*((sqr(v.x())+sqr(v.y())+sqr(v.z()))/2.5-1.0)
               -1.0/p*(s.xx()*sqr(v.x())+s.yy()*sqr(v.y())+s.zz()*sqr(v.z())+2.0*s.xy()*v.x()*v.y()+2.0*s.xz()*v.x()*v.z()+2.0*s.yz()*v.y()*v.z());

    } while(A*cloud_.rndGen().sample01<scalar>()>gamma);

    return U+u0*v;

}

Foam::vector Foam::ChapmanEnskogRelaxation::samplePostRelaxationVelocityESBGK
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

