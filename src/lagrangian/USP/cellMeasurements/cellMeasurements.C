/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "cellMeasurements.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellMeasurements::cellMeasurements
(
    const polyMesh& mesh,
    uspCloud& cloud
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud)
{}


Foam::cellMeasurements::cellMeasurements
(
    const polyMesh& mesh,
    uspCloud& cloud,
    const bool dummy
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    typeIds_(),
    binCoeff_(),
    relCoeff_(),
    rhoNMean_(),
    rhoNInstantaneous_(),
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
    nColls_(mesh_.nCells(), 0.0),
    vibrationalETotal_(),
    electronicETotal_(),
    nParcels_(),
    nParcelsXnParticle_(),
    mccSpecies_(),
    nGroundElectronicLevel_(),
    nFirstElectronicLevel_()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellMeasurements::createFields()
{

    typeIds_.setSize(cloud_.typeIdList().size());
    forAll(typeIds_, iD)
    {
        typeIds_[iD] = iD;
    }

    binCoeff_.setSize(cloud_.typeIdList().size());
    for (auto& f : binCoeff_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    relCoeff_.setSize(cloud_.typeIdList().size());
    for (auto& f : relCoeff_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    rhoNMean_.setSize(cloud_.typeIdList().size());
    for (auto& f : rhoNMean_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    rhoNInstantaneous_.setSize(cloud_.typeIdList().size());
    for (auto& f : rhoNInstantaneous_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    rhoNMeanXnParticle_.setSize(cloud_.typeIdList().size());
    for (auto& f : rhoNMeanXnParticle_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    rhoNMeanInt_.setSize(cloud_.typeIdList().size());
    for (auto& f : rhoNMeanInt_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    molsElec_.setSize(cloud_.typeIdList().size());
    for (auto& f : molsElec_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    rhoMMean_.setSize(cloud_.typeIdList().size());
    for (auto& f : rhoMMean_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    rhoMMeanXnParticle_.setSize(cloud_.typeIdList().size());
    for (auto& f : rhoMMeanXnParticle_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    linearKEMean_.setSize(cloud_.typeIdList().size());
    for (auto& f : linearKEMean_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    linearKEMeanXnParticle_.setSize(cloud_.typeIdList().size());
    for (auto& f : linearKEMeanXnParticle_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    rotationalEMean_.setSize(cloud_.typeIdList().size());
    for (auto& f : rotationalEMean_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    rotationalDofMean_.setSize(cloud_.typeIdList().size());
    for (auto& f : rotationalDofMean_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    muu_.setSize(cloud_.typeIdList().size());
    for (auto& f : muu_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    muv_.setSize(cloud_.typeIdList().size());
    for (auto& f : muv_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    muw_.setSize(cloud_.typeIdList().size());
    for (auto& f : muw_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    mvv_.setSize(cloud_.typeIdList().size());
    for (auto& f : mvv_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    mvw_.setSize(cloud_.typeIdList().size());
    for (auto& f : mvw_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    mww_.setSize(cloud_.typeIdList().size());
    for (auto& f : mww_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    mcc_.setSize(cloud_.typeIdList().size());
    for (auto& f : mcc_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    mccu_.setSize(cloud_.typeIdList().size());
    for (auto& f : mccu_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    mccv_.setSize(cloud_.typeIdList().size());
    for (auto& f : mccv_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    mccw_.setSize(cloud_.typeIdList().size());
    for (auto& f : mccw_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    eu_.setSize(cloud_.typeIdList().size());
    for (auto& f : eu_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    ev_.setSize(cloud_.typeIdList().size());
    for (auto& f : ev_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    ew_.setSize(cloud_.typeIdList().size());
    for (auto& f : ew_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    e_.setSize(cloud_.typeIdList().size());
    for (auto& f : e_)
    {
        f.setSize(mesh_.nCells(), 0.0);
    }

    momentumMean_.setSize(cloud_.typeIdList().size());
    for (auto& f : momentumMean_)
    {
        f.setSize(mesh_.nCells(), vector::zero);
    }

    momentumMeanXnParticle_.setSize(cloud_.typeIdList().size());
    for (auto& f : momentumMeanXnParticle_)
    {
        f.setSize(mesh_.nCells(), vector::zero);
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

    forAll(vibrationalETotal_, i)
    {
        vibrationalETotal_[i].setSize
        (
            cloud_.constProps(typeIds_[i]).vibrationalDoF()
        );

        forAll(vibrationalETotal_[i], j)
        {
            vibrationalETotal_[i][j].setSize(mesh_.nCells(), 0.0);
        }
    }

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


void Foam::cellMeasurements::clean()
{
    // Clean geometric fields
    forAll(typeIds_, iD)
    {

        forAll(mesh_.cells(), cell)
        {

            binCoeff_[iD][cell] = 0.0;
            relCoeff_[iD][cell] = 0.0;
            rhoNMean_[iD][cell] = 0.0;
            rhoNInstantaneous_[iD][cell] =  0.0;
            rhoNMeanXnParticle_[iD][cell] =  0.0;
            rhoNMeanInt_[iD][cell] =  0.0;
            molsElec_[iD][cell] =  0.0;
            rhoMMean_[iD][cell] =  0.0;
            rhoMMeanXnParticle_[iD][cell] =  0.0;
            linearKEMean_[iD][cell] =  0.0;
            linearKEMeanXnParticle_[iD][cell] =  0.0;
            rotationalEMean_[iD][cell] =  0.0;
            rotationalDofMean_[iD][cell] =  0.0;

            muu_[iD][cell] =  0.0;
            muv_[iD][cell] =  0.0;
            muw_[iD][cell] =  0.0;
            mvv_[iD][cell] =  0.0;
            mvw_[iD][cell] =  0.0;
            mww_[iD][cell] =  0.0;
            mcc_[iD][cell] =  0.0;
            mccu_[iD][cell] =  0.0;
            mccv_[iD][cell] =  0.0;
            mccw_[iD][cell] =  0.0;

            eu_[iD][cell] =  0.0;
            ev_[iD][cell] =  0.0;
            ew_[iD][cell] =  0.0;
            e_[iD][cell] =  0.0;

            momentumMean_[iD][cell] =  vector::zero;
            momentumMeanXnParticle_[iD][cell] = vector::zero;
            nColls_[cell] =  0.0;

            electronicETotal_[iD][cell] = 0.0;
            mccSpecies_[iD][cell] = 0.0;
            nParcels_[iD][cell] = 0.0;
            nGroundElectronicLevel_[iD][cell] = 0.0;
            nFirstElectronicLevel_[iD][cell] = 0.0;
            nParcelsXnParticle_[iD][cell] = 0.0;

            forAll(vibrationalETotal_[iD], v)
            {
                vibrationalETotal_[iD][v][cell] = 0.0;
            }

        }
    }

}


void Foam::cellMeasurements::calculateFields()
{

    const scalar deltaT = mesh_.time().deltaTValue();
    const scalar nParticle = cloud_.nParticle();

    // Calaculate USP and DSMC coefficients
    forAll(mesh_.cells(), cell)
    {
        if (cloud_.cellCollModel(cell) == cloud_.binCollModel())
        {
            forAll(typeIds_, iD)
            {
                binCoeff_[iD][cell] += deltaT;
            }       
        }
        else
        {
            forAll(typeIds_, iD)
            {
                relCoeff_[iD][cell] += deltaT;
            }
        }
    }

    // Calculate parcel properties sums
    forAllConstIters(cloud_, iter)
    {
        const uspParcel& p = iter();
        const label iD = typeIds_.find(p.typeId());

        if (iD != -1)
        {
            const label cell = p.cell();

            const scalar& CWF = p.CWF();
            const scalar& RWF = p.RWF();
            const scalar& mass = cloud_.constProps(p.typeId()).mass();
            const scalar massByMagUsq = mass*magSqr(p.U());
            const scalarList& electronicEnergies = cloud_.constProps(typeIds_[iD]).electronicEnergyList();
            const label& rotationalDof = cloud_.constProps(p.typeId()).rotationalDoF();
            const scalar& xVel = p.U().x();
            const scalar& yVel = p.U().y();
            const scalar& zVel = p.U().z();
            const vector& U = p.U();

            scalarList EVib(cloud_.constProps(typeIds_[iD]).vibrationalDoF());
            
            if (EVib.size() > 0)
            {
                forAll(EVib, i)
                {
                    EVib[i] =
                        p.vibLevel()[i]
                    *physicoChemical::k.value()
                    *cloud_.constProps(p.typeId()).thetaV()[i];

                    vibrationalETotal_[iD][i][cell] +=
                        p.vibLevel()[i]
                    *physicoChemical::k.value()
                    *cloud_.constProps(p.typeId()).thetaV()[i];
                }
            }
     
            rhoNMean_[iD][cell] += 1.0;
            rhoNInstantaneous_[iD][cell] += 1.0;
            rhoMMean_[iD][cell] += mass;
            linearKEMean_[iD][cell] += mass*(U & U);
            momentumMean_[iD][cell] += mass*U;
            rotationalEMean_[iD][cell] += p.ERot();
            rotationalDofMean_[iD][cell] += rotationalDof;
            electronicETotal_[iD][cell] +=
                electronicEnergies[p.ELevel()];
            nParcels_[iD][cell] += 1.0;
            mccSpecies_[iD][cell] += massByMagUsq;

            nParcelsXnParticle_[iD][cell] += CWF*RWF*nParticle;
            rhoNMeanXnParticle_[iD][cell] += CWF*RWF*nParticle;
            rhoMMeanXnParticle_[iD][cell] += mass*CWF*RWF*nParticle;
            momentumMeanXnParticle_[iD][cell] += mass*(U)*CWF*RWF*nParticle;
            linearKEMeanXnParticle_[iD][cell] += mass*(U & U)*CWF*RWF*nParticle;

            muu_[iD][cell] += mass*sqr(xVel);
            muv_[iD][cell] += mass*(xVel*yVel);
            muw_[iD][cell] += mass*(xVel*zVel);
            mvv_[iD][cell] += mass*sqr(yVel);
            mvw_[iD][cell] += mass*(yVel*zVel);
            mww_[iD][cell] += mass*sqr(zVel);
            mcc_[iD][cell] += massByMagUsq;
            mccu_[iD][cell] += massByMagUsq*(xVel);
            mccv_[iD][cell] += massByMagUsq*(yVel);
            mccw_[iD][cell] += massByMagUsq*(zVel);

            scalar vibEn = 0.0;

            forAll(EVib, v)
            {
                vibEn += EVib[v];
            }

            eu_[iD][cell] += (p.ERot() + vibEn)*xVel;
            ev_[iD][cell] += (p.ERot() + vibEn)*yVel;
            ew_[iD][cell] += (p.ERot() + vibEn)*zVel;
            e_[iD][cell] += p.ERot() + vibEn;

            if (rotationalDof > VSMALL)
            {
                rhoNMeanInt_[iD][cell] += 1.0;
            }

            label nElecLevels =
                cloud_.constProps(p.typeId()).nElectronicLevels();

            if (nElecLevels > 1)
            {
                molsElec_[iD][cell] += 1.0;

                if (p.ELevel() == 0)
                {
                    nGroundElectronicLevel_[iD][cell]++;
                }
                if (p.ELevel() == 1)
                {
                    nFirstElectronicLevel_[iD][cell]++;
                }
            }

        }
    }

}


// ************************************************************************* //
