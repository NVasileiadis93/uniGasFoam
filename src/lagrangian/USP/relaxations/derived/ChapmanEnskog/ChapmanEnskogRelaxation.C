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
    infoCounter_(0),
    relaxationFrequency_
    (
        IOobject
        (
            "relaxationFrequency",
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
        dimensionedVector(dimLength/dimTime, Zero),
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
        dimensionedVector(dimensionSet(1, 0, -3, 0, 0), Zero),
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
    stressTensor_
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
    typeIds_(),
    rhoNMean_(mesh_.nCells(),Zero),
    rhoNInstantaneous_(mesh_.nCells(),Zero),
    rhoNMeanXnParticle_(mesh_.nCells(),Zero),
    rhoNMeanInt_(mesh_.nCells(),Zero),
    molsElec_(mesh_.nCells(),Zero),
    rhoMMean_(mesh_.nCells(),Zero),
    rhoMMeanXnParticle_(mesh_.nCells(),Zero),
    linearKEMean_(mesh_.nCells(),Zero),
    linearKEMeanXnParticle_(mesh_.nCells(),Zero),
    rotationalEMean_(mesh_.nCells(),Zero),
    rotationalDofMean_(mesh_.nCells(),Zero),
    muu_(mesh_.nCells(),Zero),
    muv_(mesh_.nCells(),Zero),
    muw_(mesh_.nCells(),Zero),
    mvv_(mesh_.nCells(),Zero),
    mvw_(mesh_.nCells(),Zero),
    mww_(mesh_.nCells(),Zero),
    mcc_(mesh_.nCells(),Zero),
    mccu_(mesh_.nCells(),Zero),
    mccv_(mesh_.nCells(),Zero),
    mccw_(mesh_.nCells(),Zero),
    eu_(mesh_.nCells(),Zero),
    ev_(mesh_.nCells(),Zero),
    ew_(mesh_.nCells(),Zero),
    e_(mesh_.nCells(),Zero),
    totalvDof_(mesh_.nCells(),Zero),
    momentumMean_(mesh.nCells(), Zero),
    momentumMeanXnParticle_(mesh.nCells(), Zero),
    vibrationalETotal_(),
    electronicETotal_(),
    nParcels_(),
    nParcelsXnParticle_(),
    mccSpecies_(),
    vibT_(),
    vDof_(),
    nGroundElectronicLevel_(),
    nFirstElectronicLevel_()
{

    typeIds_.setSize(cloud_.typeIdList().size());
    forAll(typeIds_,iD)
    {
        typeIds_[iD] = iD;
    }

    // Outer list is typeIds, inner list is number of cells on the mesh
    vibT_.setSize(typeIds_.size());

    for (auto& v : vibT_)
    {
        v.setSize(mesh_.nCells());
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

    vDof_.setSize(typeIds_.size());

    for (auto& v : vDof_)
    {
        v.setSize(mesh_.nCells());
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

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ChapmanEnskogRelaxation::calculateProperties()
{

    scalarField vibT(mesh_.nCells(), Zero);
    scalarField vibTForOverallT(mesh_.nCells(), Zero);
    scalarField Cp(mesh_.nCells(), Zero);
    scalarField Cv(mesh_.nCells(), Zero);
    scalarField molecularMass(mesh_.nCells(), Zero);
    scalarField Cv_p(mesh_.nCells(), Zero);
    scalarField totalvDof(mesh_.nCells(), Zero);
    scalarField totalvDofOverall(mesh_.nCells(), Zero);

    const scalar nParticle = cloud_.nParticle();

    const List<DynamicList<uspParcel*>>&
        cellOccupancy = cloud_.cellOccupancy();

    forAll(cellOccupancy, cell)
    {
        const DynamicList<uspParcel*>& cellParcels(cellOccupancy[cell]);

        forAll(cellParcels, i)
        {
            const uspParcel& p = *cellParcels[i];
            const label iD = typeIds_.find(p.typeId());

            const scalar mass = cloud_.constProps(p.typeId()).mass();
            const scalar massByMagUsq = mass*magSqr(p.U());
            const scalar xVel = p.U().x();
            const scalar yVel = p.U().y();
            const scalar zVel = p.U().z();
            const vector& U = p.U();
            const scalarList& electronicEnergies =
                cloud_.constProps(typeIds_[iD]).electronicEnergyList();
            const label rotationalDof =
                cloud_.constProps(p.typeId()).rotationalDoF();
            
            scalarList EVib
            (
                cloud_.constProps(typeIds_[iD]).vibrationalDoF()
            );

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

            rhoNMean_[cell] += 1.0;
            rhoNInstantaneous_[cell] += 1.0;
            rhoMMean_[cell] += mass;
            linearKEMean_[cell] += mass*(U & U);
            momentumMean_[cell] += mass*U;
            rotationalEMean_[cell] += p.ERot();
            rotationalDofMean_[cell] += rotationalDof;
            electronicETotal_[iD][cell] +=
                electronicEnergies[p.ELevel()];
            nParcels_[iD][cell] += 1.0;
            mccSpecies_[iD][cell] += massByMagUsq;

            scalar CWF = 
                cloud_.cellWF(p.cell());
            scalar RWF =
                cloud_.axiRWF(p.position());

            nParcelsXnParticle_[iD][cell] += CWF*RWF*nParticle;
            rhoNMeanXnParticle_[cell] += CWF*RWF*nParticle;
            rhoMMeanXnParticle_[cell] += mass*CWF*RWF*nParticle;
            momentumMeanXnParticle_[cell] += mass*(U)*CWF*RWF*nParticle;
            linearKEMeanXnParticle_[cell] += mass*(U & U)*CWF*RWF*nParticle;

            muu_[cell] += mass*sqr(xVel);
            muv_[cell] += mass*(xVel*yVel);
            muw_[cell] += mass*(xVel*zVel);
            mvv_[cell] += mass*sqr(yVel);
            mvw_[cell] += mass*(yVel*zVel);
            mww_[cell] += mass*sqr(zVel);
            mcc_[cell] += massByMagUsq;
            mccu_[cell] += massByMagUsq*(xVel);
            mccv_[cell] += massByMagUsq*(yVel);
            mccw_[cell] += massByMagUsq*(zVel);

            scalar vibEn = 0.0;

            forAll(EVib, v)
            {
                vibEn += EVib[v];
            }

            eu_[cell] += (p.ERot() + vibEn)*xVel;
            ev_[cell] += (p.ERot() + vibEn)*yVel;
            ew_[cell] += (p.ERot() + vibEn)*zVel;
            e_[cell] += p.ERot() + vibEn;

            if (rotationalDof > VSMALL)
            {
                rhoNMeanInt_[cell] += 1.0;
            }

            label nElecLevels =
                cloud_.constProps(p.typeId()).nElectronicLevels();

            if (nElecLevels > 1)
            {
                molsElec_[cell] += 1.0;

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

    forAll(cellOccupancy, cell)
    {

        if (rhoNMean_[cell] > VSMALL)
        {
            const scalar cellVolume = mesh_.cellVolumes()[cell];

            rhoN_[cell] =
                (rhoNMeanXnParticle_[cell])/(cellVolume);

            rhoM_[cell] =
                (rhoMMeanXnParticle_[cell])/(cellVolume);

            scalar rhoMMean =
                rhoMMeanXnParticle_[cell]/(cellVolume);
            UMean_[cell] =
                momentumMeanXnParticle_[cell]
               /(rhoMMean*cellVolume);
            scalar linearKEMean =
                0.5*linearKEMeanXnParticle_[cell]
               /(cellVolume);
            scalar rhoNMean =
                rhoNMeanXnParticle_[cell]/(cellVolume);

            translationalT_[cell] =
                2.0/(3.0*physicoChemical::k.value()*rhoNMean)
               *(
                    linearKEMean
                  - 0.5*rhoMMean*(UMean_[cell] & UMean_[cell])
                );

            p_[cell] =
                rhoN_[cell]*physicoChemical::k.value()
                *translationalT_[cell];

                if (rotationalDofMean_[cell] > VSMALL)
                {
                    scalar rotationalEMean =
                        rotationalEMean_[cell];
                    scalar rotationalDofMean =
                        rotationalDofMean_[cell];

                    rotationalT_[cell] =
                        (2.0/physicoChemical::k.value())
                       *(rotationalEMean/rotationalDofMean);
                }
                else
                {
                    rotationalT_[cell] = 0.0;
                }

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

                scalar nRotDof = 0.0;

                if (rhoNMean_[cell] > VSMALL)
                {
                    nRotDof = rotationalDofMean_[cell] / rhoNMean_[cell];
                }

                overallT_[cell] =
                    (
                        (3.0*translationalT_[cell])
                      + (nRotDof*rotationalT_[cell])
                      + (totalvDof_[cell]*vibrationalT_[cell])
                      + (totalEDof*electronicT_[cell])
                    )
                   /(3.0 + nRotDof + totalvDof_[cell] + totalEDof);

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

                stressTensor_[cell] = -pressureTensor_[cell];
                stressTensor_[cell].xx() += scalarPressure;
                stressTensor_[cell].yy() += scalarPressure;
                stressTensor_[cell].zz() += scalarPressure;

                //terms involving pressure tensor should not be
                //multiplied by the number density (see Bird corrigendum)

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
            UMean_[cell] = vector::zero;
            translationalT_[cell] = 0.0;
            p_[cell] = 0.0;
            pressureTensor_[cell] = tensor::zero;
            stressTensor_[cell] = tensor::zero;
            heatFluxVector_[cell] = vector::zero;
        }

    }

    //collision frequency
    scalar Tref = cloud_.constProps(0).Tref();
    forAll(cellOccupancy, cell)
    {
        const DynamicList<uspParcel*>& cellParcels(cellOccupancy[cell]);

        scalar parcels = 0.0;
        scalar omegaTotal = 0.0;
        scalar totalViscRef = 0.0;

        forAll(cellParcels, i)
        {

            uspParcel& p = *cellParcels[i];
            const scalar mass = cloud_.constProps(p.typeId()).mass();
            const scalar omega = cloud_.constProps(p.typeId()).omega();
            const scalar d = cloud_.constProps(p.typeId()).d();

            parcels++;
            omegaTotal += omega;
            totalViscRef += 7.5*sqrt(mass*constant::physicoChemical::k.value()*Tref)
                            /(sqrt(constant::mathematical::pi)*(5.0-2.0*omega)*(7.0-2.0*omega)*sqr(d));

        }

        scalar omegaMean = omegaTotal/parcels;
        scalar viscRefMean = totalViscRef/parcels;
        relaxationFrequency_[cell] = Pr_*p_[cell]/(viscRefMean*pow(translationalT_[cell]/Tref,omegaMean));

    }


    forAll(cellOccupancy, cell)
    {
        std::cout << "--------------rhoM-rhoM-P-T-Pr-visc-----------------" << std::endl;
        std::cout << rhoM_[cell] << " " << rhoN_[cell] << " " << p_[cell] << " " << translationalT_[cell]  << std::endl;
        std::cout << "------------------------U--------------------------" << std::endl;
        std::cout << UMean_[cell].x() << " " << UMean_[cell].y() << " " << UMean_[cell].z() << std::endl;
        std::cout << "------------------------pT-------------------------" << std::endl;
        std::cout << pressureTensor_[cell].xx() << " " << pressureTensor_[cell].xy() << " " << pressureTensor_[cell].xz() << std::endl;
        std::cout << pressureTensor_[cell].yx() << " " << pressureTensor_[cell].yy() << " " << pressureTensor_[cell].yz() << std::endl;
        std::cout << pressureTensor_[cell].zx() << " " << pressureTensor_[cell].zy() << " " << pressureTensor_[cell].zz() << std::endl;
        std::cin.get();
    }

    // reset 
    forAll(rhoNMean_, cell)
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
    }

    forAll(electronicETotal_, iD)
    {
        forAll(electronicETotal_[iD], cell)
        {
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

void Foam::ChapmanEnskogRelaxation::relax()
{

    const scalar deltaT = cloud_.mesh().time().deltaTValue();

    // Calculate required macroscopic properties
    calculateProperties();

    // Relax particles to local distribution
    label relaxations = 0;

    List<DynamicList<uspParcel*>>&
        cellOccupancy = cloud_.cellOccupancy();

    forAll(cellOccupancy, cell)
    {
        const DynamicList<uspParcel*>& cellParcels(cellOccupancy[cell]);

        forAll(cellParcels, i)
        {

            if (cloud_.rndGen().sample01<scalar>() < 1.0-exp(-relaxationFrequency_[cell]*deltaT)) //cloud_.rndGen().sample01<scalar>()
            {

                uspParcel& p = *cellParcels[i];
                const scalar mass = cloud_.constProps(p.typeId()).mass();

                // Calculate breakdown parameter
                scalar mostProbableSpeed(
                    cloud_.maxwellianMostProbableSpeed
                    (
                        translationalT_[cell],
                        mass
                    ));

                scalar maxHeatFlux = -1.0;
                forAll(heatFluxVector_[cell],i) 
                {
                    if (maxHeatFlux < fabs(heatFluxVector_[cell][i])) 
                    {
                        maxHeatFlux = fabs(heatFluxVector_[cell][i]);
                    }
                }
                maxHeatFlux = 2.0*maxHeatFlux/(p_[cell]*mostProbableSpeed);

                scalar maxStress = -1.0;
                forAll(stressTensor_[cell],i) 
                {
                    if (maxStress < fabs(stressTensor_[cell][i])) 
                    {
                        maxStress = fabs(stressTensor_[cell][i]);
                    }
                }
                maxStress = maxStress/p_[cell];

                const scalar& breakdownParameter = max(maxHeatFlux,maxStress);

                // Relax particle
                p.U() = samplePostRelaxationVelocity
                        (
                            breakdownParameter,
                            mass,
                            p_[cell],
                            translationalT_[cell],
                            UMean_[cell],
                            heatFluxVector_[cell],
                            stressTensor_[cell]
                        );

                relaxations++;
            
            }
        }

        const scalar cellVolume = mesh_.cellVolumes()[cell];
        scalar nParticle = cloud_.nParticle();
        scalar newTranslationalT;
        vector newUMean;
        scalar rhoNMean = 0.0;
        scalar rhoMMean = 0.0;
        scalar linearKEMean = 0.0;
        scalar rhoNMeanXnParticle = 0.0;
        scalar rhoMMeanXnParticle = 0.0;
        scalar linearKEMeanXnParticle = 0.0;
        vector momentumMean = vector::zero;
        vector momentumMeanXnParticle = vector::zero;

        forAll(cellParcels, i)
        {

            const uspParcel& p = *cellParcels[i];

            const scalar& mass = cloud_.constProps(p.typeId()).mass();
            const vector& U = p.U();

            rhoNMean += 1.0;
            rhoMMean += mass;
            linearKEMean += mass*(U & U);
            momentumMean += mass*U;

            scalar CWF = 
                cloud_.cellWF(p.cell());
            scalar RWF =
                cloud_.axiRWF(p.position());

            rhoNMeanXnParticle += CWF*RWF*nParticle;
            rhoMMeanXnParticle += mass*CWF*RWF*nParticle;
            momentumMeanXnParticle += mass*(U)*CWF*RWF*nParticle;
            linearKEMeanXnParticle += mass*(U & U)*CWF*RWF*nParticle;
        }

        newUMean = momentumMean/rhoMMean;
        linearKEMean = 0.5*linearKEMeanXnParticle/(cellVolume);
        newTranslationalT =
            2.0/(3.0*physicoChemical::k.value()*rhoNMeanXnParticle/(cellVolume))
           *(
                linearKEMean
              - 0.5*rhoMMeanXnParticle/(cellVolume)*(newUMean & newUMean)
            );

        forAll(cellParcels, i)
        {
                uspParcel& p = *cellParcels[i];
                p.U() = UMean_[cell] + (p.U()-newUMean)*sqrt(translationalT_[cell]/newTranslationalT);
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

    vector pU;
    scalar gamma;
    scalar A = 1.0 + 60.0*B;
    scalar u0(cloud_.maxwellianMostProbableSpeed(T,m));
    
    do {

        pU=cloud_.rndGen().GaussNormal<vector>()/sqrt(2.0);

        gamma=1.0+(2.0/u0*(q.x()*pU.x()+q.y()*pU.y()+q.z()*pU.z())*(2.0/5.0*(pow(pU.x(),2)+pow(pU.y(),2)+pow(pU.z(),2))-1.0)
             -2.0*(s.xy()*pU.x()*pU.y()+s.xz()*pU.x()*pU.z()+s.yz()*pU.y()*pU.z())
             -s.xx()*(pow(pU.x(),2)-pow(pU.z(),2))-s.yy()*(pow(pU.y(),2)-pow(pU.z(),2)))/p;

    } while(A*cloud_.rndGen().sample01<scalar>()>gamma);

    return U+u0*pU;

}

// ************************************************************************* //

