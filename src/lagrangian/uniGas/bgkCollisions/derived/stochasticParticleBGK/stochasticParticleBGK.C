/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "stochasticParticleBGK.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(stochasticParticleBGK, 0);

addToRunTimeSelectionTable(bgkCollisionModel, stochasticParticleBGK, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stochasticParticleBGK::stochasticParticleBGK
(
    const dictionary& dict,
    const polyMesh& mesh,
    uniGasCloud& cloud
)
:
    bgkCollisionModel(dict, mesh, cloud),
    propertiesDict_(dict.subDict("collisionProperties")),
    Tref_(propertiesDict_.get<scalar>("Tref")),
    macroInterpolation_(propertiesDict_.getOrDefault<bool>("macroInterpolation", false)),
    infoCounter_(0),
    shufflePasses_(5),
    performCollision_(mesh_.nCells(), true),
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
    collFreq_
    (
        IOobject
        (
            "collFreq",
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

        forAll(vibrationalETotal_[i], j)
        {
            vibrationalETotal_[i][j].setSize(mesh_.nCells(), 0.0);
        }
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::stochasticParticleBGK::calculateProperties()
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

            // Scale macroscopic properties
            if (rhoNMean_[cell] > 1.0)
            {
                p_[cell] = rhoNMean_[cell]/(rhoNMean_[cell]-1.0)*p_[cell];
                translationalT_[cell] = rhoNMean_[cell]/(rhoNMean_[cell]-1.0)*translationalT_[cell];
            }
            else
            {
                performCollision_[cell] = false;
            }

        }
        else
        {
            performCollision_[cell] = false;
            rhoN_[cell] = 0.0;
            p_[cell] = 0.0;
            translationalT_[cell] = 0.0;
            UMean_[cell] = vector::zero;
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
        scalar vibT = 0.0;
        scalar totalvDof = 0.0;
        scalarList degreesOfFreedomSpecies(typeIds_.size(), 0.0);
        scalarList vibTID(vibrationalETotal_.size(), 0.0);
        List<scalarList> dofMode;
        List<scalarList> vibTMode;
        dofMode.setSize(typeIds_.size());
        vibTMode.setSize(typeIds_.size());

        forAll(dofMode, iD)
        {
            const auto& constProp = cloud_.constProps(typeIds_[iD]);
            
            dofMode[iD].setSize(constProp.vibrationalDoF(), 0.0);
            vibTMode[iD].setSize(constProp.vibrationalDoF(), 0.0);
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
                    const auto& constProp = 
                        cloud_.constProps(typeIds_[iD]);
                    
                    scalar thetaV = constProp.thetaV()[v];

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

            totalvDof += degreesOfFreedomSpecies[iD];

            if
            (
                rhoNMeanInt_[cell] > VSMALL
             && rhoNMean_[cell] > VSMALL
             && nParcels_[iD][cell] > VSMALL
            )
            {
                vibT += vibTID[iD]*nParcels_[iD][cell]/rhoNMeanInt_[cell];;
            }
        }

        vibrationalT_[cell] = vibT;

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
              + (totalvDof*vibrationalT_[cell])
              + (totalEDof*electronicT_[cell])
            )
           /(3.0 + nRotDof + totalvDof + totalEDof);

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

                Prandtl_[cell] += nParcels_[iD][cell]*(5+rotDoF)/(7.5+rotDoF);

            }
            viscosity /= rhoNMean_[cell];
            Prandtl_[cell] /= rhoNMean_[cell]; 

            collFreq_[cell] = p_[cell]/viscosity;
        }
        else
        {
            performCollision_[cell] = false;
            Prandtl_[cell] = 0.0;
            collFreq_[cell] = 0.0;
        }

    }

    //Correct boundary conditions
    Prandtl_.correctBoundaryConditions();
    collFreq_.correctBoundaryConditions();
    rhoN_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    translationalT_.correctBoundaryConditions();
    rotationalT_.correctBoundaryConditions();
    vibrationalT_.correctBoundaryConditions();
    electronicT_.correctBoundaryConditions();
    overallT_.correctBoundaryConditions();
    UMean_.correctBoundaryConditions();

}

void Foam::stochasticParticleBGK::resetProperties()
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

        performCollision_[cell] = true;

    }

}

void Foam::stochasticParticleBGK::collide()
{

    const scalar& deltaT = cloud_.mesh().time().deltaTValue();

    // Collide particles to local distribution
    label collisions = 0;

    List<DynamicList<uniGasParcel*>>& cellOccupancy = cloud_.cellOccupancy();

    // Calculate required macroscopic properties
    calculateProperties();

    // Create macroscopic quantities interpolations
    autoPtr <Foam::interpolation<scalar>> translationalTInterp;
    autoPtr <Foam::interpolation<vector>> UMeanInterp;

    if (macroInterpolation_)
    {
        translationalTInterp = Foam::interpolationCellPoint<scalar>::New("cellPoint", translationalT_);
        UMeanInterp = Foam::interpolationCellPoint<vector>::New("cellPoint", UMean_);
    }

    forAll(cellOccupancy, cell)
    {

        if (cloud_.cellCollModelId()[cell] == cloud_.bgkCollModelId() && performCollision_[cell])
        {

            DynamicList<uniGasParcel*>& cellParcels(cellOccupancy[cell]);

            for (label pass = 1; pass <= shufflePasses_; pass++)
            {
                cloud_.rndGen().shuffle(cellParcels);
            }

            scalar parcelCollisions = rhoNMean_[cell]*(1.0-exp(-collFreq_[cell]*deltaT));
            label nParcelCollisions = int(parcelCollisions);
            if (cloud_.rndGen().sample01<scalar>() < (parcelCollisions-nParcelCollisions))
            {
                nParcelCollisions++;
            }
            nParcelCollisions = min(nParcelCollisions,rhoNMean_[cell]);

            for(label i = 0; i < nParcelCollisions; i++)
            {

                uniGasParcel& parcel = *cellParcels[i];
                const scalar& mass = cloud_.constProps(parcel.typeId()).mass();
                const vector& position = parcel.position();

                // Collide particle
                if (macroInterpolation_)
                {
                    parcel.U() = samplePostCollisionVelocity
                            (
                                mass,
                                translationalTInterp().interpolate(position, cell),
                                UMeanInterp().interpolate(position, cell)
                            );
                }
                else
                {
                    parcel.U() = samplePostCollisionVelocity
                            (
                                mass,
                                translationalT_[cell],
                                UMean_[cell]
                            );                    
                }

                collisions++;
                
            }

            conserveMomentumAndEnergy(cell);

        }
    }

    resetProperties();

    infoCounter_++;
    reduce(collisions, sumOp<label>());

    if (infoCounter_ >= cloud_.nTerminalOutputs())
    {
        if (collisions>0)
        {
            Info<< "    BGK collisions                  = "
                << collisions << nl
                << endl;
            infoCounter_ = 0;
        }
        else
        {
            Info<< "    No BGK collisions" << endl;
            infoCounter_ = 0;
        }
    }

}

void Foam::stochasticParticleBGK::conserveMomentumAndEnergy
(
    const label& cell
)
{

    // Energy conservation scheme
    const scalar nParticle = cloud_.nParticle();
    //const scalar cellVolume = mesh_.cellVolumes()[cell];
    const DynamicList<uniGasParcel*>& cellParcels(cloud_.cellOccupancy()[cell]);

    scalar linearKEMeanXnParticle = 0.0;
    vector momentumMeanXnParticle = vector::zero;

    forAll(cellParcels, i)
    {

        const uniGasParcel& p = *cellParcels[i];

        const scalar mass = cloud_.constProps(p.typeId()).mass();
        const vector& U = p.U();
        const scalar& CWF = p.CWF();
        const scalar& RWF = p.RWF();

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
                uniGasParcel& p = *cellParcels[i];
                p.U() = UMean_[cell] + (p.U()-postUMean)*sqrt(translationalT_[cell]/postTranslationalT);
        }
    }

}

Foam::vector Foam::stochasticParticleBGK::samplePostCollisionVelocity
(   
    const scalar& m,
    const scalar& T,
    const vector& U
)
{

    scalar u0(cloud_.maxwellianMostProbableSpeed(T,m));

    vector v = cloud_.rndGen().GaussNormal<vector>()/sqrt(2.0);

    return U+u0*v;

}

const Foam::dictionary&
Foam::stochasticParticleBGK::propertiesDict() const
{
    return propertiesDict_;
}

// ************************************************************************* //

