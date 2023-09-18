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

#include "uspZone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspZone, 0);
addToRunTimeSelectionTable(uspField, uspZone, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspZone::uspZone
(
    const Time& t,
    const polyMesh& mesh,
    uspCloud& cloud,
    const dictionary& dict
)
:
    uspField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    regionName_(propsDict_.get<word>("zone")),
    regionId_(-1),
    sampleInterval_(propsDict_.getOrDefault("sampleInterval", 1)),
    sampleCounter_(0),
    fieldName_(propsDict_.get<word>("field")),
    totalVolume_(0.0),
    typeIds_(),
    timeIndex_(0),
    averagingCounter_(0.0),

    mols_(0.0),
    uspMols_(0.0),
    molsInt_(0.0),
    molsElec_(0.0),
    mass_(0.0),
    mcc_(0.0),
    mom_(Zero),
    UCollected_(Zero),
    rotationalEMean_(0.0),
    rotationalDofMean_(0.0),

    muu_(Zero),
    mccu_(Zero),
    nColls_(0),

    eu_(Zero),
    e_(0.0),
    stepIndex_(0),
    speciesMols_(),
    mccSpecies_(),
    electronicETotal_(),
    nParticlesGroundElectronicState_(),
    nParticlesFirstElectronicState_(),
    vDof_(),
    mfp_(),
    mcr_(),

    N_(),
    rhoN_(),
    rhoM_(),
    UMean_(),
    translationalTemperature_(),
    rotationalTemperature_(),
    vibrationalTemperature_(),
    electronicTemperature_(),
    overallTemperature_(),
    scalarPressure_(),
    pField_(),
    tauField_(),
    qField_(),
    qInternalField_(),
    qTranslationalField_(),
    meanFreePath_(),
    meanCollisionRate_(),
    meanCollisionTime_(),
    meanCollisionTimeTimeStepRatio_(),
    Ma_(),
    vibrationalETotal_(),

    outputField_(4, true),
    instantaneous_(propsDict_.getOrDefault<bool>("instantaneous", false)),
    averagingAcrossManyRuns_
    (
        propsDict_.getOrDefault<bool>("averagingAcrossManyRuns", false)
    ),
    measureCollisionRateOnly_
    (
        propsDict_.getOrDefault<bool>("measureCollisionRateOnly", false)
    )
{
    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if (regionId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find region: " << regionName_ << nl << "in: "
            << mesh_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

    typeIds_ = cloud_.getTypeIDs(propsDict_);

    // Set the total volume
    const labelList& cells = cellZones[regionId_];

    for (label celli : cells)
    {
        totalVolume_ += mesh_.cellVolumes()[celli];
    }

    if (Pstream::parRun())
    {
        reduce(totalVolume_, sumOp<scalar>());
    }

    // ---------------------------------------------------

    // Instantaneous
    const scalar deltaT = timeVel_.mdTimeInterval().deltaT();
    scalar writeInterval = t.controlDict().get<scalar>("writeInterval");
    label nBins = label(writeInterval/deltaT);
    nSteps_ = 1;

    if (!instantaneous_) // cumulative
    {
        nBins = 1;
        nSteps_ = label(writeInterval/deltaT);
    }

    N_.setSize(nBins, 0.0);
    rhoN_.setSize(nBins, 0.0);
    rhoM_.setSize(nBins, 0.0);
    UMean_.setSize(nBins, Zero);
    UCAM_.setSize(nBins, Zero);
    translationalTemperature_.setSize(nBins, 0.0);
    rotationalTemperature_.setSize(nBins, 0.0);
    vibrationalTemperature_.setSize(nBins, 0.0);
    electronicTemperature_.setSize(nBins, 0.0);
    overallTemperature_.setSize(nBins, 0.0);
    scalarPressure_.setSize(nBins, 0.0);
    pField_.setSize(nBins, Zero);
    tauField_.setSize(nBins, Zero);
    qField_.setSize(nBins, Zero);
    qInternalField_.setSize(nBins, Zero);
    qTranslationalField_.setSize(nBins, Zero);
    vibrationalModeTemperatures_.setSize(nBins, Zero);
    meanFreePath_.setSize(nBins, 0.0);
    meanCollisionRate_.setSize(nBins, 0.0);
    meanCollisionTime_.setSize(nBins, 0.0);
    meanCollisionTimeTimeStepRatio_.setSize(nBins, 0.0);
    measuredCollisionRate_.setSize(nBins, 0.0);
    Ma_.setSize(nBins, 0.0);

    speciesMols_.setSize(typeIds_.size(), 0.0);
    mccSpecies_.setSize(typeIds_.size(), 0.0);
    vibrationalETotal_.setSize(typeIds_.size());
    electronicETotal_.setSize(typeIds_.size(), 0.0);
    nParticlesGroundElectronicState_.setSize(typeIds_.size(), 0.0);
    nParticlesFirstElectronicState_.setSize(typeIds_.size(), 0.0);
    vDof_.setSize(typeIds_.size(), 0.0);
    mfp_.setSize(typeIds_.size(), 0.0);
    mcr_.setSize(typeIds_.size(), 0.0);

    // Read stored data from dictionary
    if (averagingAcrossManyRuns_)
    {
        Info<< nl << "Averaging across many runs initiated." << nl << endl;

        read();
    }

    // Read in stored data from dictionary
    if (measureCollisionRateOnly_)
    {
        Info<< nl << "Only measuring and outputting"
            << " collision rate. " << nl << endl;

        read();
    }

    // Choice of measurement property to output

    if (propsDict_.found("outputProperties"))
    {
        const wordHashSet measurements(propsDict_.lookup("outputProperties"));

        outputField_[0] = measurements.found("density");
        outputField_[1] = measurements.found("velocity");
        outputField_[2] = measurements.found("temperature");
        outputField_[3] = measurements.found("pressure");

        // Check for user:
        forAllConstIters(measurements, iter)
        {
            const word& propertyName(iter());

            if
            (
                (propertyName != "density")
             && (propertyName != "velocity")
             && (propertyName != "temperature")
             && (propertyName != "pressure")
            )
            {
                FatalErrorInFunction
                    << "Cannot find measurement property: " << propertyName
                    << nl << "in: "
                    << mesh_.time().system()/"fieldPropertiesDict"
                    << exit(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspZone::read()
{
    localIOdictionary dict
    (
        IOobject
        (
            "zone_" + fieldName_ + "_" + regionName_,
            mesh_.time().timeName(),
            "uniform",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    dict.readIfPresent("mols", mols_);
    dict.readIfPresent("uspMols", uspMols_);
    dict.readIfPresent("molsInt", molsInt_);
    dict.readIfPresent("molsElec", molsElec_);
    dict.readIfPresent("mass", mass_);
    dict.readIfPresent("mcc", mcc_);
    dict.readIfPresent("mccSpecies", mccSpecies_);
    dict.readIfPresent("mom", mom_);
    dict.readIfPresent("UCollected", UCollected_);
    dict.readIfPresent("rotationalEMean", rotationalEMean_);
    dict.readIfPresent("rotationalDofMean", rotationalDofMean_);

    dict.readIfPresent("muu", muu_);

    dict.readIfPresent("mccu", mccu_);
    dict.readIfPresent("nColls", nColls_);
    dict.readIfPresent("eu", eu_);
    dict.readIfPresent("e", e_);

    dict.readIfPresent("vibrationalETotal", vibrationalETotal_);
    dict.readIfPresent("electronicETotal", electronicETotal_);
    dict.readIfPresent
    (
        "nParticlesGroundElectronicState",
        nParticlesGroundElectronicState_
    );
    dict.readIfPresent
    (
        "nParticlesFirstElectronicState",
        nParticlesFirstElectronicState_
    );
    dict.readIfPresent("speciesMols", speciesMols_);

    dict.readIfPresent("averagingCounter", averagingCounter_);
}


void Foam::uspZone::write()
{
    if (mesh_.time().writeTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "zone_" + fieldName_ + "_" + regionName_,
                mesh_.time().timeName(),
                "uniform",
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );

        dict.add("mols", mols_);
        dict.add("uspMols", uspMols_);
        dict.add("molsInt", molsInt_);
        dict.add("molsElec", molsElec_);
        dict.add("mass", mass_);
        dict.add("mcc", mcc_);
        dict.add("mccSpecies", mccSpecies_);
        dict.add("mom", mom_);
        dict.add("UCollected", UCollected_);
        dict.add("rotationalEMean", rotationalEMean_);
        dict.add("rotationalDofMean", rotationalDofMean_);

        dict.add("muu", muu_);

        dict.add("mccu", mccu_);
        dict.add("nColls", nColls_);
        dict.add("eu", eu_);
        dict.add("e", e_);

        dict.add("vibrationalETotal", vibrationalETotal_);
        dict.add("electronicETotal", electronicETotal_);
        dict.add
        (
            "nParticlesGroundElectronicState",
            nParticlesGroundElectronicState_
        );
        dict.add
        (
            "nParticlesFirstElectronicState",
            nParticlesFirstElectronicState_
        );
        dict.add("speciesMols", speciesMols_);

        dict.add("averagingCounter", averagingCounter_);

        dict.regIOobject::writeObject
        (
            IOstreamOption(mesh_.time().writeFormat()),
            true
        );
    }
}


void Foam::uspZone::createField()
{
    forAll(vibrationalETotal_, iD)
    {
        vibrationalETotal_[iD].setSize
        (
            cloud_.constProps(typeIds_[iD]).vibrationalDoF(),
            scalar(0)
        );
    }
}


void Foam::uspZone::calculateField()
{
    ++sampleCounter_;

    if (sampleInterval_ <= sampleCounter_)
    {
        averagingCounter_ += 1.0;

        const auto& cellOccupancy = cloud_.cellOccupancy();

        const labelList& cells = mesh_.cellZones()[regionId_];

        for (const label celli : cells)
        {
            for (const uspParcel* p : cellOccupancy[celli])
            {
                label iD = typeIds_.find(p->typeId());

                if (iD != -1)
                {

                    scalar CWF = cloud_.cellWF(p->cell());
                    scalar RWF = cloud_.axiRWF(p->position());
                    scalar nParticle = cloud_.nParticle()*CWF*RWF;
                    
                    nColls_ += cloud_.cellPropMeasurements().nColls()[celli]*nParticle;
                    
                    const auto& constProp = cloud_.constProps(p->typeId());

                    const label rotationalDof = constProp.rotationalDoF();

                    const scalar mass = constProp.mass()*nParticle;
                    const scalarList& electronicEnergies =
                        constProp.electronicEnergyList();

                    mols_ += nParticle;
                    uspMols_ += 1.0;
                    mass_ += mass;
                    mcc_ += mass*magSqr(p->U());
                    UCollected_ += p->U();
                    mom_ += mass*p->U();
                    rotationalEMean_ += p->ERot();
                    rotationalDofMean_ += rotationalDof;

                    muu_ += mass*symm(p->U()*p->U());
                    mccu_ += mass*magSqr(p->U())*p->U();

                    scalarList EVib
                    (
                        constProp.vibrationalDoF(), 0.0
                    );

                    forAll(EVib, i)
                    {
                        EVib[i] =
                            p->vibLevel()[i]*physicoChemical::k.value()
                           *constProp.thetaV()[i];
                    }

                    scalar vibE = sum(EVib);

                    eu_ += nParticle*(p->ERot() + vibE)*p->U();
                    e_ += nParticle*(p->ERot() + vibE);

                    vibrationalETotal_[iD] += EVib;
                    electronicETotal_[iD] += electronicEnergies[p->ELevel()];
                    speciesMols_[iD] += 1.0;
                    mccSpecies_[iD] += mass*magSqr(p->U());

                    if (rotationalDof > VSMALL)
                    {
                        molsInt_ += 1.0;
                    }

                    if (constProp.nElectronicLevels() > 1)
                    {
                        molsElec_ += 1.0;
                    }

                    if (p->ELevel() == 0)
                    {
                        nParticlesGroundElectronicState_[iD] += 1.0;
                    }

                    if (p->ELevel() == 1)
                    {
                        nParticlesFirstElectronicState_[iD] += 1.0;
                    }
                }
            }
        }
        sampleCounter_ = 0;
    }

    ++stepIndex_;

    if (stepIndex_ >= nSteps_)
    {
        stepIndex_ = 0;

        scalar mass = mass_;
        scalar mols = mols_;
        scalar uspMols = uspMols_;
        scalar molsInt = molsInt_;
        scalar molsElec = molsElec_;
        scalar mcc = mcc_;
        vector mom = mom_;
        vector UCollected = UCollected_;
        scalar rotationalEMean = rotationalEMean_;
        scalar rotationalDofMean = rotationalDofMean_;

        symmTensor muu = muu_;
        vector mccu = mccu_;
        scalar nColls = nColls_;

        vector eu = eu_;
        scalar e = e_;

        List<scalarField> vibrationalETotal(vibrationalETotal_);
        scalarField electronicETotal(electronicETotal_);
        scalarField nParticlesGroundElectronicState
        (
            nParticlesGroundElectronicState_
        );
        scalarField nParticlesFirstElectronicState
        (
            nParticlesFirstElectronicState_
        );
        scalarField speciesMols(speciesMols_);
        scalarField mccSpecies(mccSpecies_);

        // parallel communication

        if (Pstream::parRun())
        {
            reduce(mols, sumOp<scalar>());
            reduce(uspMols, sumOp<scalar>());
            reduce(molsInt, sumOp<scalar>());
            reduce(molsElec, sumOp<scalar>());
            reduce(mass, sumOp<scalar>());
            reduce(mcc, sumOp<scalar>());
            reduce(mom, sumOp<vector>());
            reduce(UCollected, sumOp<vector>());
            reduce(rotationalEMean, sumOp<scalar>());
            reduce(rotationalDofMean, sumOp<scalar>());

            reduce(muu, sumOp<symmTensor>());
            reduce(mccu, sumOp<vector>());
            reduce(nColls, sumOp<scalar>());

            reduce(eu, sumOp<vector>());
            reduce(e, sumOp<scalar>());

            forAll(vibrationalETotal, iD)
            {
                reduce(vibrationalETotal[iD], sumOp<scalarList>());
                reduce(electronicETotal[iD], sumOp<scalar>());
                reduce(speciesMols[iD], sumOp<scalar>());
                reduce(mccSpecies[iD], sumOp<scalar>());
                reduce(nParticlesGroundElectronicState[iD], sumOp<scalar>());
                reduce(nParticlesFirstElectronicState[iD], sumOp<scalar>());
            }
        }

        const scalar volume = totalVolume_;
        label n = timeIndex_;

        N_[n] = uspMols/averagingCounter_;
        rhoN_[n] = mols/(averagingCounter_*volume);
        rhoM_[n] = mass/(averagingCounter_*volume);

        if (uspMols > 0)
        {
            UMean_[n] = UCollected/uspMols;

            UCAM_[n] = mom/mass;

            translationalTemperature_[n] =
                (1.0/(3.0*physicoChemical::k.value()))
               *(mcc/mols - mass/mols*magSqr(UMean_[n]));


            if (rotationalDofMean > VSMALL)
            {
                rotationalTemperature_[n] =
                    (2.0/physicoChemical::k.value())
                   *(rotationalEMean/rotationalDofMean);
            }
            else
            {
                rotationalTemperature_[n] = 0.0;
            }

            pField_[n] = rhoN_[n]*(muu/mols - mass/mols*symm(sqr(UMean_[n])));

            scalarPressure_[n] = (1.0/3.0)*tr(pField_[n]);

            // Make reference
            tauField_[n] = -pField_[n] + scalarPressure_[n]*I;

            // Terms involving pressure tensor should not be multiplied
            // by the number density (see Bird corrigendum)
            qField_[n] =
                rhoN_[n]
               *(
                    0.5*mccu/mols
                  - 0.5*mcc/mols*UMean_[n]
                  + eu/mols
                  - e/mols*UMean_[n]
                )
              - (pField_[n] & UMean_[n]);

            qInternalField_[n] = rhoN_[n]*(eu/mols - e/mols*UMean_[n]);

            qTranslationalField_[n] =
                rhoN_[n]*(0.5*mccu/mols - 0.5*mcc/mols*UMean_[n])
              - (pField_[n] & UMean_[n]);

            // Vibrational temperature
            scalar totalvDof = 0.0;
            scalar totalvDofOverall = 0.0;
            scalar vibT = 0.0;
            scalarList degreesOfFreedomSpecies(typeIds_.size(), 0.0);
            scalarList vibTID(vibrationalETotal.size(),0.0);

            List<scalarList> degreesOfFreedomMode(typeIds_.size());
            List<scalarList> vibTMode(typeIds_.size());

            forAll(degreesOfFreedomMode, iD)
            {
                degreesOfFreedomMode[iD].setSize
                (
                    vibrationalETotal.size(),
                    Zero
                );
                vibTMode[iD].setSize(vibrationalETotal.size(), Zero);
            }

            forAll(vibrationalETotal, iD)
            {
                forAll(vibrationalETotal[iD], m)
                {
                    if
                    (
                        vibrationalETotal[iD][m] > VSMALL
                     && speciesMols[iD] > VSMALL
                    )
                    {
                        const auto& constProp = cloud_.constProps(typeIds_[iD]);
                        
                        scalar thetaV = constProp.thetaV()[m];

                        scalarList vibrationalEMean =
                            vibrationalETotal[iD]/speciesMols[iD];

                        scalar iMean = 0.0;

                        iMean =
                            vibrationalEMean[m]
                           /(physicoChemical::k.value()*thetaV);

                        vibTMode[iD][m] = thetaV/log(1.0 + (1.0/iMean));

                        degreesOfFreedomMode[iD][m] =
                            (2.0*thetaV/vibTMode[iD][m])
                           /(exp(thetaV/vibTMode[iD][m]) - 1.0);
                    }
                }

                forAll(vibrationalETotal[iD], m)
                {
                    degreesOfFreedomSpecies[iD] += degreesOfFreedomMode[iD][m];
                }

                forAll(vibrationalETotal[iD], m)
                {
                    if (degreesOfFreedomSpecies[iD] > VSMALL)
                    {
                        vibTID[iD] +=
                            vibTMode[iD][m]
                           *degreesOfFreedomMode[iD][m]
                           /degreesOfFreedomSpecies[iD];
                    }
                }

                totalvDof += degreesOfFreedomSpecies[iD];

                scalar fraction = 0.0;

                if (molsInt > VSMALL)
                {
                    fraction = speciesMols[iD]/molsInt;
                }

                scalar fractionOverall = speciesMols[iD]/mols;

                if (fraction > SMALL)
                {
                    totalvDofOverall += totalvDof*(fractionOverall/fraction);
                }

                vibT += vibTID[iD]*fraction;
            }

            vibrationalTemperature_[n] = vibT;
            vibrationalModeTemperatures_[n].x() = vibTMode[0][0];
            vibrationalModeTemperatures_[n].y() = vibTMode[0][1];
            vibrationalModeTemperatures_[n].z() = vibTMode[0][2];

            // Electronic temperature
            scalar totalEDof = 0.0;
            scalar elecT = 0.0;

            forAll(speciesMols, iD)
            {
                const auto& cpi = cloud_.constProps(typeIds_[iD]);

                label nElectronicLevels = cpi.nElectronicLevels();

                if
                (
                    nElectronicLevels > 1
                 && speciesMols[iD] > VSMALL
                 && molsElec > VSMALL
                )
                {
                    const scalarList& electronicEnergies =
                        cpi.electronicEnergyList();
                    const labelList& degeneracies = cpi.degeneracyList();

                    if
                    (
                        nParticlesGroundElectronicState[iD] > VSMALL
                     && nParticlesFirstElectronicState[iD] > VSMALL
                     && (
                            nParticlesGroundElectronicState[iD]*degeneracies[1]
                         != nParticlesFirstElectronicState[iD]*degeneracies[0]
                        )
                    )
                    {
                        scalar fraction = speciesMols[iD]/molsElec;

                        scalar elecTID =
                            (electronicEnergies[1] - electronicEnergies[0])
                           /(
                                physicoChemical::k.value()
                               *log
                                (
                                    (
                                        nParticlesGroundElectronicState[iD]
                                       *degeneracies[1]
                                    )
                                   /(
                                       nParticlesFirstElectronicState[iD]
                                      *degeneracies[0]
                                    )
                                )
                            );

                        if (elecTID > VSMALL)
                        {
                            elecT += fraction*elecTID;
                        }

                        scalar eDof =
                            2*electronicETotal[iD]/speciesMols[iD]
                           /(physicoChemical::k.value()*elecTID);

                        totalEDof += fraction*eDof;
                    }
                }
            }

            electronicTemperature_[n] = elecT;

            // Overall temperature

            scalar nRotDof = 0.0;

            if (uspMols > VSMALL)
            {
                nRotDof = rotationalDofMean/uspMols;
            }

            overallTemperature_[n] =
                (
                    3*translationalTemperature_[n]
                  + nRotDof*rotationalTemperature_[n]
                  + totalvDofOverall*vibrationalTemperature_[n]
                  + totalEDof*electronicTemperature_[n]
                )
               /(3.0 + nRotDof + totalvDofOverall + totalEDof);

            forAll(mfp_, iD)
            {
                label qspec = 0;

                const auto& cpi = cloud_.constProps(typeIds_[iD]);

                for (qspec=0; qspec<typeIds_.size(); ++qspec)
                {                   
                    const auto& cpq = cloud_.constProps(typeIds_[qspec]);

                    scalar dPQ = 0.5*(cpi.d() + cpq.d());

                    scalar omegaPQ = 0.5*(cpi.omega() + cpq.omega());

                    scalar massRatio = cpi.mass()/cpq.mass();

                    scalar reducedMass =
                        (cpi.mass()*cpq.mass())/(cpi.mass() + cpq.mass());

                    if
                    (
                        speciesMols[qspec] > VSMALL
                     && translationalTemperature_[n] > VSMALL
                    )
                    {
                        scalar nDensQ =
                            (cloud_.nParticle()*speciesMols[qspec])
                           /(volume*averagingCounter_);

                        // Bird, eq (4.76)
                        mfp_[iD] +=
                            pi*dPQ*dPQ*nDensQ
                           *pow
                            (
                                273.0/translationalTemperature_[n],
                                omegaPQ - 0.5
                            )
                            *sqrt(1.0 + massRatio);

                        // Bird, eq (4.74)
                        mcr_[iD] +=
                            2.0*sqrt(pi)*dPQ*dPQ*nDensQ
                           *pow
                            (
                                translationalTemperature_[n]/273.0,
                                1.0 - omegaPQ
                            )
                           *sqrt
                            (
                                2.0*physicoChemical::k.value()*273.0/reducedMass
                            );
                    }
                }

                if (mfp_[iD] > VSMALL)
                {
                    mfp_[iD] = 1.0/mfp_[iD];
                }
            }

            meanFreePath_[n] = 0.0;
            meanCollisionRate_[n] = 0.0;
            meanCollisionTime_[n] = 0.0;
            meanCollisionTimeTimeStepRatio_[n] = 0.0;

            if (rhoN_[n] > VSMALL)
            {
                forAll(mfp_, iD)
                {
                    scalar nDensP =
                        (cloud_.nParticle()*speciesMols[iD])
                       /(volume*averagingCounter_);

                    // Bird, eq (4.77)
                    meanFreePath_[n] += mfp_[iD]*nDensP/rhoN_[n];

                    // Bird, eq (1.38)
                    meanCollisionRate_[n] += mcr_[iD]*nDensP/rhoN_[n];
                }
            }
            else
            {
                meanFreePath_[n] = GREAT;
                meanCollisionRate_[n] = 0.0;
            }

            const scalar deltaT = mesh_.time().deltaTValue();

            if (meanCollisionRate_[n] > VSMALL)
            {
                meanCollisionTime_[n] = 1.0/meanCollisionRate_[n];
                meanCollisionTimeTimeStepRatio_[n] =
                    meanCollisionTime_[n]/deltaT;
            }
            else
            {
                meanCollisionTime_[n] = GREAT;
                meanCollisionTimeTimeStepRatio_[n] = GREAT;
            }

            measuredCollisionRate_[n] =
                nColls/(volume*deltaT*averagingCounter_);

            mfp_ = scalar(0);
            mcr_ = scalar(0);

            scalar molecularMass = 0.0;
            scalar molarconstantPressureSpecificHeat = 0.0;
            scalar molarconstantVolumeSpecificHeat = 0.0;
            scalar speedOfSound = 0.0;
            scalar gasConstant = 0.0;
            scalar gamma = 0.0;

            forAll(mfp_, iD)
            {
                if (uspMols > VSMALL)
                {
                    const auto& cpi = cloud_.constProps(typeIds_[iD]);
                    molecularMass += cpi.mass()*(speciesMols[iD]/uspMols);
                    molarconstantPressureSpecificHeat +=
                        (5.0 + cpi.rotationalDoF())
                       *(speciesMols[iD]/uspMols);
                    molarconstantVolumeSpecificHeat +=
                        (3.0 + cpi.rotationalDoF())
                       *(speciesMols[iD]/uspMols);
                }
            }

            if (molecularMass > VSMALL)
            {
                // R = k/m
                gasConstant = physicoChemical::k.value()/molecularMass;
            }

            if (molarconstantVolumeSpecificHeat > VSMALL)
            {
                // gamma = cP/cV
                gamma =
                    molarconstantPressureSpecificHeat
                   /molarconstantVolumeSpecificHeat;
            }

            if
            (
                translationalTemperature_[n] > VSMALL
             && gamma > VSMALL
             && gasConstant > VSMALL
            )
            {
                speedOfSound =
                    sqrt(gamma*gasConstant*translationalTemperature_[n]);
            }

            if (speedOfSound > VSMALL)
            {
                Ma_[n] = mag(UMean_[n])/speedOfSound;
            }
            else
            {
                Ma_[n] = 0.0;
            }
        }

        if (timeVel_.resetFieldsAtOutput() && (mesh_.time().value() < timeVel_.resetFieldsAtOutputUntilTime()+0.5*mesh_.time().deltaTValue()))
        {
            // Reset fields
            averagingCounter_ = 0.0;

            mols_ = 0.0;
            uspMols_ = 0.0;
            molsInt_ = 0.0;
            molsElec_ = 0.0;
            mass_ = 0.0;
            mcc_ = 0.0;
            mom_ = Zero;
            UCollected_ = Zero;
            rotationalEMean_ = 0.0;
            rotationalDofMean_ = 0.0;

            muu_ = Zero;
            mccu_ = Zero;
            mcc_ = 0.0;
            nColls_ = 0.0;
            eu_ = Zero;
            e_ = 0.0;
            speciesMols_ = 0.0;
            mccSpecies_ = 0.0;

            forAll(electronicETotal_, iD)
            {
                electronicETotal_[iD] = 0.0;
                nParticlesGroundElectronicState_[iD] = 0.0;
                nParticlesFirstElectronicState_[iD] = 0.0;

                forAll(vibrationalETotal_[iD], j)
                {
                    vibrationalETotal_[iD][j] = 0.0;
                }
            }
        }

        if (averagingAcrossManyRuns_)
        {
            write();
        }

        ++timeIndex_;
    }
}


void Foam::uspZone::writeField()
{
    const Time& runTime = mesh_.time();

    if (runTime.writeTime())
    {
        timeIndex_ = 0;

        if (Pstream::master())
        {
            scalarField timeField(N_.size(), Zero);
            const scalar deltaT = timeVel_.mdTimeInterval().deltaT();

            forAll(timeField, t)
            {
                timeField[N_.size()-t-1] = runTime.timeOutputValue() - deltaT*t;
            }

            const word prefix("zone_" + regionName_ + "_" + fieldName_);

            if (measureCollisionRateOnly_)
            {
                writeTimeData
                (
                    casePath_,
                    prefix + "_measuredCollisionRate.xyz",
                    timeField,
                    measuredCollisionRate_,
                    true
                );
            }
            else
            {
                // Output densities
                if (outputField_[0])
                {
                    writeTimeData
                    (
                        casePath_,
                        prefix + "_N.xy",
                        timeField,
                        N_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_rhoN.xy",
                        timeField,
                        rhoN_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_rhoM.xy",
                        timeField,
                        rhoM_,
                        true
                    );
                }

                // output velocities
                if (outputField_[1])
                {
                    writeTimeData
                    (
                        casePath_,
                        prefix + "_U_SAM.xyz",
                        timeField,
                        UMean_,
                        true
                    );


                    writeTimeData
                    (
                        casePath_,
                        prefix + "_U_CAM.xyz",
                        timeField,
                        UCAM_,
                        true
                    );

                }

                // output temperature
                if (outputField_[2])
                {
                    writeTimeData
                    (
                        casePath_,
                        prefix + "_translationalTemperature.xy",
                        timeField,
                        translationalTemperature_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_rotationalTemperature.xy",
                        timeField,
                        rotationalTemperature_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_vibrationalTemperature.xy",
                        timeField,
                        vibrationalTemperature_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_electronicTemperature.xy",
                        timeField,
                        electronicTemperature_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_overallTemperature.xy",
                        timeField,
                        overallTemperature_,
                        true
                    );
                }

                // output pressure
                if (outputField_[3])
                {
                    writeTimeData
                    (
                        casePath_,
                        prefix + "_pressureTensor.xyz",
                        timeField,
                        pField_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_p.xy",
                        timeField,
                        scalarPressure_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_stressTensor.xyz",
                        timeField,
                        tauField_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_heatFluxVector.xyz",
                        timeField,
                        qField_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_internalHeatFluxVector.xyz",
                        timeField,
                        qInternalField_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_translationalHeatFluxVector.xyz",
                        timeField,
                        qTranslationalField_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_variableHardSphereMeanFreePath.xyz",
                        timeField,
                        meanFreePath_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_meanCollisionRate.xyz",
                        timeField,
                        meanCollisionRate_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_meanCollisionTime.xyz",
                        timeField,
                        meanCollisionTime_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_meanCollisionTimeTimeStepRatio.xyz",
                        timeField,
                        meanCollisionTimeTimeStepRatio_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_measuredCollisionRate.xyz",
                        timeField,
                        measuredCollisionRate_,
                        true
                    );

                    writeTimeData
                    (
                        casePath_,
                        prefix + "_Ma.xyz",
                        timeField,
                        Ma_,
                        true
                    );
                }
            }
        }
    }
}


void Foam::uspZone::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uspField::updateProperties(dict);
}


// ************************************************************************* //