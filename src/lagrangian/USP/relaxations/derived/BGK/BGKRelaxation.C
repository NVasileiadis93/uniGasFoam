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

#include "BGKRelaxation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(BGKRelaxation, 0);

addToRunTimeSelectionTable(relaxationModel, BGKRelaxation, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BGKRelaxation::BGKRelaxation
(
    const dictionary& dict,
    const polyMesh& mesh,
    uspCloud& cloud
)
:
    relaxationModel(dict, mesh, cloud),
    infoCounter_(0),
    steps_(0),
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
    typeIds_(),
    rhoNMean_(mesh_.nCells(),Zero),
    rhoNMeanXnParticle_(mesh_.nCells(),Zero),
    rhoMMean_(mesh_.nCells(),Zero),
    rhoMMeanXnParticle_(mesh_.nCells(),Zero),
    linearKEMean_(mesh_.nCells(),Zero),
    linearKEMeanXnParticle_(mesh_.nCells(),Zero),
    momentumMean_(mesh.nCells(), Zero),
    momentumMeanXnParticle_(mesh.nCells(), Zero)
{

    typeIds_.setSize(cloud_.typeIdList().size());
    forAll(typeIds_,iD)
    {
        typeIds_[iD] = iD;
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::BGKRelaxation::calculateProperties()
{

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
            const vector& U = p.U();

            rhoNMean_[cell] += 1.0;
            rhoMMean_[cell] += mass;
            linearKEMean_[cell] += mass*(U & U);
            momentumMean_[cell] += mass*U;

            scalar CWF = 
                cloud_.cellWF(p.cell());
            scalar RWF =
                cloud_.axiRWF(p.position());

            rhoNMeanXnParticle_[cell] += CWF*RWF*nParticle;
            rhoMMeanXnParticle_[cell] += mass*CWF*RWF*nParticle;
            momentumMeanXnParticle_[cell] += mass*(U)*CWF*RWF*nParticle;
            linearKEMeanXnParticle_[cell] += mass*(U & U)*CWF*RWF*nParticle;

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

        }
        else
        {
            rhoN_[cell] = 0.0;
            rhoM_[cell] = 0.0;
            UMean_[cell] = vector::zero;
            translationalT_[cell] = 0.0;
            p_[cell] = 0.0;
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
        relaxationFrequency_[cell] = p_[cell]/(viscRefMean*pow(translationalT_[cell]/Tref,omegaMean));

    }

    // reset 
    forAll(rhoNMean_, cell)
    {
        rhoNMean_[cell] = scalar(0.0);
        rhoMMean_[cell] = scalar(0.0);
        linearKEMean_[cell] = scalar(0.0);
        momentumMean_[cell] = vector::zero;
        rhoNMeanXnParticle_[cell] = scalar(0.0);
        rhoMMeanXnParticle_[cell] = scalar(0.0);
        momentumMeanXnParticle_[cell] = vector::zero;
        linearKEMeanXnParticle_[cell] = scalar(0.0);
    }

}

void Foam::BGKRelaxation::relax()
{

    steps_++;

    // Calculate required macroscopic properties
    calculateProperties();

    // Relax particles to local distribution
    label relaxations = 0;

    const scalar deltaT = cloud_.mesh().time().deltaTValue();

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

                // Relax particle
                p.U() = samplePostRelaxationVelocity
                        (
                            mass,
                            translationalT_[cell],
                            UMean_[cell]
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

            const scalar mass = cloud_.constProps(p.typeId()).mass();
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


Foam::vector Foam::BGKRelaxation::samplePostRelaxationVelocity
(   
    scalar m,
    scalar T,
    vector U
)
{

    return cloud_.maxwellianMostProbableSpeed(T,m)/sqrt(2.0)*cloud_.rndGen().GaussNormal<vector>() + U;

}

// ************************************************************************* //

