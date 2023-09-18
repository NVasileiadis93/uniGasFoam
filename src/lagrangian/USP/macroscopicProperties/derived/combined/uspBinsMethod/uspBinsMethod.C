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

#include "uspBinsMethod.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspBinsMethod, 0);
addToRunTimeSelectionTable(uspField, uspBinsMethod, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspBinsMethod::uspBinsMethod
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
    fieldName_(propsDict_.get<word>("field")),
    typeIds_(),
    sampleInterval_(1),
    sampleCounter_(0),
    averagingCounter_(0.0),

    volumeOfCellsInBin_(),
    mols_(),
    uspMols_(),
    molsInt_(),
    molsElec_(),
    mass_(),
    mcc_(),
    UCollected_(),
    rotationalEMean_(),
    rotationalDofMean_(),

    muu_(),
    muv_(),
    muw_(),
    mvv_(),
    mvw_(),
    mww_(),
    scalarPressure_(),
    mccu_(),
    mccv_(),
    mccw_(),

    eu_(),
    ev_(),
    ew_(),
    e_(),

    speciesMols_(),
    mccSpecies_(),
    electronicETotal_(),
    nParticlesGroundElectronicState_(),
    nParticlesFirstElectronicState_(),
    vDof_(),
    mfp_(),

    binVolume_(),
    N_(),
    M_(),
    rhoN_(),
    rhoM_(),
    averageMass_(),
    UMean_(),
    UCAM_(),
    scalarUMean_(),
    translationalTemperature_(),
    rotationalTemperature_(),
    vibrationalTemperature_(),
    electronicTemperature_(),
    overallTemperature_(),

    pField_(),
    tauField_(),
    qField_(),
    qInternalField_(),
    qTranslationalField_(),
    meanFreePath_(),
    Ma_(),
    vibrationalETotal_(),

    outputField_(4, true),
    averagingAcrossManyRuns_(false),
    permeabilityMeasurements_(false),
    rotationalMeasurements_(false),
    wallVelocity_(0.0),
    wallRadius_(0.0)
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

    // Create bin model
    binModel_ = binModel::New(mesh, propsDict_);

    const label nBins = binModel_->nBins();

    volumeOfCellsInBin_.setSize(nBins, Zero);
    mols_.setSize(nBins, Zero);
    uspMols_.setSize(nBins, Zero);
    molsInt_.setSize(nBins, Zero);
    molsElec_.setSize(nBins, Zero);
    mass_.setSize(nBins, Zero);
    mcc_.setSize(nBins, Zero);
    mom_.setSize(nBins, Zero);
    UCollected_.setSize(nBins, Zero);
    rotationalEMean_.setSize(nBins, Zero);
    rotationalDofMean_.setSize(nBins, Zero);

    muu_.setSize(nBins, Zero);
    muv_.setSize(nBins, Zero);
    muw_.setSize(nBins, Zero);
    mvv_.setSize(nBins, Zero);
    mvw_.setSize(nBins, Zero);
    mww_.setSize(nBins, Zero);

    mccu_.setSize(nBins, Zero);
    mccv_.setSize(nBins, Zero);
    mccw_.setSize(nBins, Zero);

    eu_.setSize(nBins, Zero);
    ev_.setSize(nBins, Zero);
    ew_.setSize(nBins, Zero);
    e_.setSize(nBins, Zero);

    vibrationalETotal_.setSize(nBins);
    electronicETotal_.setSize(nBins);
    nParticlesGroundElectronicState_.setSize(nBins);
    nParticlesFirstElectronicState_.setSize(nBins);
    speciesMols_.setSize(nBins);
    mccSpecies_.setSize(nBins);
    vDof_.setSize(nBins);
    mfp_.setSize(nBins);

    binVolume_.setSize(nBins, Zero);
    N_.setSize(nBins, Zero);
    M_.setSize(nBins, Zero);
    rhoN_.setSize(nBins, Zero);
    rhoM_.setSize(nBins, Zero);
    averageMass_.setSize(nBins, Zero);
    UMean_.setSize(nBins, Zero);
    UCAM_.setSize(nBins, Zero);
    scalarUMean_.setSize(nBins, Zero);
    translationalTemperature_.setSize(nBins, Zero);
    rotationalTemperature_.setSize(nBins, Zero);
    vibrationalTemperature_.setSize(nBins, Zero);
    electronicTemperature_.setSize(nBins, Zero);
    overallTemperature_.setSize(nBins, Zero);

    scalarPressure_.setSize(nBins, Zero);
    pField_.setSize(nBins, Zero);
    tauField_.setSize(nBins, Zero);
    qField_.setSize(nBins, Zero);
    qInternalField_.setSize(nBins, Zero);
    qTranslationalField_.setSize(nBins, Zero);
    meanFreePath_.setSize(nBins, Zero);
    Ma_.setSize(nBins, Zero);

    // Outer list is the bins, inner list is type ids
    forAll(electronicETotal_, b)
    {
        vibrationalETotal_[b].setSize(typeIds_.size());
        electronicETotal_[b].setSize(typeIds_.size(), Zero);
        nParticlesGroundElectronicState_[b].setSize(typeIds_.size(), Zero);
        nParticlesFirstElectronicState_[b].setSize(typeIds_.size(), Zero);
        speciesMols_[b].setSize(typeIds_.size(), Zero);
        mccSpecies_[b].setSize(typeIds_.size(), Zero);
        vDof_[b].setSize(typeIds_.size(), Zero);
        mfp_[b].setSize(typeIds_.size(), Zero);
    }

    propsDict_.readIfPresent<label>("sampleInterval", sampleInterval_);

    if
    (
        propsDict_.readIfPresent<bool>
        (
            "averagingAcrossManyRuns",
            averagingAcrossManyRuns_
        )
    )
    {
        // Read stored data from dictionary
        if (averagingAcrossManyRuns_)
        {
            Info << nl << "Averaging across many runs initiated." << nl << endl;

            readIn();
        }
    }

    if
    (
        propsDict_.readIfPresent<bool>
        (
            "permeabilityMeasurements",
            permeabilityMeasurements_
        )
    )
    {
        // Read stored data from dictionary
        if (permeabilityMeasurements_)
        {
            Info<< nl << "Measuring volume of cells in bins; "
                << "ensure your bins are aligned with the mesh!" << nl << endl;
        }
    }

    if
    (
        propsDict_.readIfPresent<bool>
        (
            "rotationalMeasurements",
            rotationalMeasurements_
        )
    )
    {
        // Read stored data from dictionary
        if (rotationalMeasurements_)
        {
            wallVelocity_ = propsDict_.get<scalar>("wallVelocity");

            wallRadius_ = propsDict_.get<scalar>("wallRadius");

            if (wallRadius_ < VSMALL)
            {
                FatalErrorInFunction
                    << "Wall radius is: " << wallRadius_
                    << ", must be greater than zero."
                    << exit(FatalError);
            }
        }
    }

    // Choice of measurement property to output

    if (propsDict_.found("outputProperties"))
    {
        const List<word> measurements(propsDict_.lookup("outputProperties"));

        DynamicList<word> propertyNames(0);

        forAll(measurements, i)
        {
            const word& propertyName(measurements[i]);

            if (propertyNames.find(propertyName) == -1)
            {
                propertyNames.append(propertyName);
            }
        }

        propertyNames.shrink();

        if (propertyNames.find("density") == -1)
        {
            outputField_[0] = false;
        }

        if (propertyNames.find("velocity") == -1)
        {
            outputField_[1] = false;
        }

        if (propertyNames.find("temperature") == -1)
        {
            outputField_[2] = false;
        }

        if (propertyNames.find("pressure") == -1)
        {
            outputField_[3] = false;
        }

        // Check for user:
        for (const word& propertyName : propertyNames)
        {
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

void Foam::uspBinsMethod::readIn()
{
    localIOdictionary dict
    (
        IOobject
        (
            "binsMethod_"+fieldName_+"_"+regionName_,
            mesh_.time().timeName(),
            "uniform",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    dict.readIfPresent("volumeOfCellsInBin", volumeOfCellsInBin_);
    dict.readIfPresent("mols", mols_);
    dict.readIfPresent("uspMols", uspMols_);
    dict.readIfPresent("molsInt", molsInt_);
    dict.readIfPresent("molsElec", molsElec_);
    dict.readIfPresent("mass", mass_);
    dict.readIfPresent("mcc", mcc_);
    dict.readIfPresent("mom", mom_);
    dict.readIfPresent("UCollected", UCollected_);
    dict.readIfPresent("rotationalEMean", rotationalEMean_);
    dict.readIfPresent("rotationalDofMean", rotationalDofMean_);

    dict.readIfPresent("muu", muu_);
    dict.readIfPresent("muv", muv_);
    dict.readIfPresent("muw", muw_);
    dict.readIfPresent("mvv", mvv_);
    dict.readIfPresent("mvw", mvw_);
    dict.readIfPresent("mww", mww_);

    dict.readIfPresent("mccu", mccu_);
    dict.readIfPresent("mccv", mccv_);
    dict.readIfPresent("mccw", mccw_);
    dict.readIfPresent("eu", eu_);
    dict.readIfPresent("ev", ev_);
    dict.readIfPresent("ew", ew_);
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
    dict.readIfPresent("mccSpecies", mccSpecies_);

    dict.readIfPresent("averagingCounter", averagingCounter_);
}


void Foam::uspBinsMethod::writeOut()
{
    if (mesh_.time().writeTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "binsMethod_"+fieldName_+"_"+regionName_,
                mesh_.time().timeName(),
                "uniform",
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );

        dict.add("volumeOfCellsInBin", volumeOfCellsInBin_);
        dict.add("mols", mols_);
        dict.add("uspMols", uspMols_);
        dict.add("molsInt", molsInt_);
        dict.add("molsElec", molsElec_);
        dict.add("mass", mass_);
        dict.add("mcc", mcc_);
        dict.add("mom", mom_);
        dict.add("UCollected", UCollected_);
        dict.add("rotationalEMean", rotationalEMean_);
        dict.add("rotationalDofMean", rotationalDofMean_);

        dict.add("muu", muu_);
        dict.add("muv", muv_);
        dict.add("muw", muw_);
        dict.add("mvv", mvv_);
        dict.add("mvw", mvw_);
        dict.add("mww", mww_);

        dict.add("mccu", mccu_);
        dict.add("mccv", mccv_);
        dict.add("mccw", mccw_);
        dict.add("eu", eu_);
        dict.add("ev", ev_);
        dict.add("ew", ew_);
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
        dict.add("mccSpecies", mccSpecies_);

        dict.add("averagingCounter", averagingCounter_);

        dict.regIOobject::writeObject
        (
            IOstreamOption(mesh_.time().writeFormat()),
            true
        );
    }
}


void Foam::uspBinsMethod::createField()
{
    for (auto& vet : vibrationalETotal_)
    {
        forAll(vet, iD)
        {
            vet[iD].setSize
            (
                cloud_.constProps(typeIds_[iD]).vibrationalDoF(),
                Zero
            );
        }
    }
}


void Foam::uspBinsMethod::calculateField()
{
    ++sampleCounter_;

    if (sampleInterval_ <= sampleCounter_)
    {
        averagingCounter_ += 1.0;

        const auto& cellOccupancy = cloud_.cellOccupancy();

        const labelList& cells = mesh_.cellZones()[regionId_];

        for (const label cellI : cells)
        {
            if (averagingCounter_ < SMALL+1.0)
            {
                const vector& cC = mesh_.cellCentres()[cellI];

                label bin = binModel_->isPointWithinBin(cC, cellI);

                if (bin != -1)
                {
                    volumeOfCellsInBin_[bin] += mesh_.cellVolumes()[cellI];
                }
            }

            const List<uspParcel*>& molsInCell = cellOccupancy[cellI];

            for (uspParcel* p : molsInCell)
            {
                label iD = typeIds_.find(p->typeId());

                const vector rI(p->position());

                label n = binModel_->isPointWithinBin(rI, cellI);

                if (n != -1 && iD != -1)
                {
                    const auto& constProp = cloud_.constProps(p->typeId());

                    scalar nParticle = cloud_.nParticle();

                    scalar CWF = cloud_.cellWF(cellI);
                    if (cloud_.axisymmetric())
                    {
                        scalar RWF = cloud_.axiRWF(p->position());
                        nParticle *= CWF*RWF;
                    }

                    const label rotationalDof = constProp.rotationalDoF();

                    const scalar mass = constProp.mass()*nParticle;
                    const scalarList& electronicEnergies =
                        cloud_.constProps(typeIds_[iD]).electronicEnergyList();

                    mols_[n] += nParticle;
                    uspMols_[n] += 1.0;
                    mass_[n] += mass;
                    mcc_[n] += mass*magSqr(p->U());
                    UCollected_[n] += p->U();
                    mom_[n] += mass*p->U();
                    rotationalEMean_[n] += p->ERot();
                    rotationalDofMean_[n] += rotationalDof;

                    muu_[n] += mass*sqr(p->U().x());
                    muv_[n] += mass*((p->U().x())*(p->U().y()));
                    muw_[n] += mass*((p->U().x())*(p->U().z()));
                    mvv_[n] += mass*sqr(p->U().y());
                    mvw_[n] += mass*((p->U().y())*(p->U().z()));
                    mww_[n] += mass*sqr(p->U().z());
                    mccu_[n] += mass*magSqr(p->U())*(p->U().x());
                    mccv_[n] += mass*magSqr(p->U())*(p->U().y());
                    mccw_[n] += mass*magSqr(p->U())*(p->U().z());


                    scalarList EVib
                    (
                        cloud_.constProps(typeIds_[iD]).vibrationalDoF()
                    );

                    forAll(EVib, i)
                    {
                        EVib[i] =
                            p->vibLevel()[i]*physicoChemical::k.value()
                           *constProp.thetaV()[i];
                    }

                    scalar vibE = sum(EVib);

                    eu_[n] += nParticle*(p->ERot() + vibE)*(p->U().x());
                    ev_[n] += nParticle*(p->ERot() + vibE)*(p->U().y());
                    ew_[n] += nParticle*(p->ERot() + vibE)*(p->U().z());
                    e_[n] += nParticle*(p->ERot() + vibE);

                    vibrationalETotal_[n][iD] += EVib;
                    electronicETotal_[n][iD] += electronicEnergies[p->ELevel()];
                    speciesMols_[n][iD] += 1.0;
                    mccSpecies_[n][iD] += mass*magSqr(p->U());

                    if (rotationalDof > VSMALL)
                    {
                        molsInt_[n] += 1.0;
                    }

                    if (constProp.nElectronicLevels() > 1)
                    {
                        molsElec_[n] += 1.0;
                    }

                    if (p->ELevel() == 0)
                    {
                        nParticlesGroundElectronicState_[n][iD] += 1.0;
                    }

                    if (p->ELevel() == 1)
                    {
                        nParticlesFirstElectronicState_[n][iD] += 1.0;
                    }
                }
            }
        }
        sampleCounter_ = 0;
    }

    if (timeVel_.averagingTime())
    {
        scalarField volumeOfCellsInBin(volumeOfCellsInBin_);
        scalarField mass(mass_);
        scalarField molsInt(molsInt_);
        scalarField molsElec(molsElec_);
        scalarField mols(mols_);
        scalarField uspMols(uspMols_);
        scalarField mcc(mcc_);
        vectorField mom(mom_);
        vectorField UCollected(UCollected_);
        scalarField rotationalEMean(rotationalEMean_);
        scalarField rotationalDofMean(rotationalDofMean_);

        scalarField muu(muu_);
        scalarField muv(muv_);
        scalarField muw(muw_);
        scalarField mvv(mvv_);
        scalarField mvw(mvw_);
        scalarField mww(mww_);
        scalarField mccu(mccu_);
        scalarField mccv(mccv_);
        scalarField mccw(mccw_);

        scalarField eu(eu_);
        scalarField ev(ev_);
        scalarField ew(ew_);
        scalarField e(e_);

        List<List<scalarField>> vibrationalETotal = vibrationalETotal_;
        List<scalarField> electronicETotal = electronicETotal_;
        List<scalarField> nParticlesGroundElectronicState =
            nParticlesGroundElectronicState_;
        List<scalarField> nParticlesFirstElectronicState =
            nParticlesFirstElectronicState_;
        List<scalarField> speciesMols = speciesMols_;
        List<scalarField> mccSpecies = mccSpecies_;

        //  parallel communication

        if (Pstream::parRun())
        {
            forAll(mols, n)
            {
                reduce(volumeOfCellsInBin[n], sumOp<scalar>());
                reduce(mols[n], sumOp<scalar>());
                reduce(uspMols[n], sumOp<scalar>());
                reduce(molsInt[n], sumOp<scalar>());
                reduce(molsElec[n], sumOp<scalar>());
                reduce(mass[n], sumOp<scalar>());
                reduce(mcc[n], sumOp<scalar>());
                reduce(mom[n], sumOp<vector>());
                reduce(UCollected[n], sumOp<vector>());
                reduce(rotationalEMean[n], sumOp<scalar>());
                reduce(rotationalDofMean[n], sumOp<scalar>());

                reduce(muu[n], sumOp<scalar>());
                reduce(muv[n], sumOp<scalar>());
                reduce(muw[n], sumOp<scalar>());
                reduce(mvv[n], sumOp<scalar>());
                reduce(mvw[n], sumOp<scalar>());
                reduce(mww[n], sumOp<scalar>());
                reduce(mccu[n], sumOp<scalar>());
                reduce(mccv[n], sumOp<scalar>());
                reduce(mccw[n], sumOp<scalar>());

                reduce(eu[n], sumOp<scalar>());
                reduce(ev[n], sumOp<scalar>());
                reduce(ew[n], sumOp<scalar>());
                reduce(e[n], sumOp<scalar>());
            }

            forAll(electronicETotal, b)
            {
                forAll(electronicETotal[b], n)
                {
                    reduce(electronicETotal[b][n], sumOp<scalar>());
                    reduce
                    (
                        nParticlesGroundElectronicState[b][n],
                        sumOp<scalar>()
                    );
                    reduce
                    (
                        nParticlesFirstElectronicState[b][n],
                        sumOp<scalar>()
                    );
                    reduce(speciesMols[b][n], sumOp<scalar>());
                    reduce(mccSpecies[b][n], sumOp<scalar>());
                }
            }

            forAll(vibrationalETotal, b)
            {
                forAll(vibrationalETotal[b], n)
                {
                    reduce(vibrationalETotal[b][n], sumOp<scalarList>());
                }
            }
        }

        forAll(mols, n)
        {
            const scalar volume = binModel_->binVolume(n);

            binVolume_[n] = volumeOfCellsInBin[n];
            N_[n] = uspMols[n]/averagingCounter_;
            M_[n] = mass[n]/averagingCounter_;
            rhoN_[n] = (mols[n])/(averagingCounter_*volume);
            rhoM_[n] = mass[n]/(averagingCounter_*volume);

            if (uspMols[n] > 0)
            {
                averageMass_[n] = M_[n]/(N_[n]*cloud_.nParticle());

                UMean_[n] = UCollected[n]/uspMols[n];

                if(rotationalMeasurements_)
                {
                    scalarUMean_[n] =
                        wallVelocity_*binModel_->binPositions()[n]/wallRadius_;
                }

                UCAM_[n] = mom[n]/mass[n];

                translationalTemperature_[n] =
                    (1.0/(3.0*physicoChemical::k.value()))
                   *(
                        mcc[n]/mols[n]
                      - mass[n]/mols[n]*magSqr(UMean_[n])
                    );

                if (rotationalMeasurements_)
                {
                    translationalTemperature_[n] =
                        (1.0/(3.0*physicoChemical::k.value()))
                       *(
                            mcc[n]/mols[n]
                          - mass[n]/mols[n]*sqr(scalarUMean_[n])
                        );
                }

                if (rotationalDofMean[n] > VSMALL)
                {
                    rotationalTemperature_[n] =
                        (2.0/physicoChemical::k.value())
                       *(rotationalEMean[n]/rotationalDofMean[n]);
                }
                else
                {
                    rotationalTemperature_[n] = 0.0;
                }

                tensorField p(uspMols.size(), Zero);

                p[n].xx() = rhoN_[n]*(muu[n]/mols[n] - (mass[n]/mols[n])*UMean_[n].x()*UMean_[n].x());
                p[n].xy() = rhoN_[n]*(muv[n]/mols[n] - (mass[n]/mols[n])*UMean_[n].x()*UMean_[n].y());
                p[n].xz() = rhoN_[n]*(muw[n]/mols[n] - (mass[n]/mols[n])*UMean_[n].x()*UMean_[n].z());
                p[n].yx() = p[n].xy();
                p[n].yy() = rhoN_[n]*(mvv[n]/mols[n] - (mass[n]/mols[n])*UMean_[n].y()*UMean_[n].y());
                p[n].yz() = rhoN_[n]*(mvw[n]/mols[n] - (mass[n]/mols[n])*UMean_[n].y()*UMean_[n].z());
                p[n].zx() = p[n].xz();
                p[n].zy() = p[n].yz();
                p[n].zz() = rhoN_[n]*(mww[n]/mols[n] - (mass[n]/mols[n])*UMean_[n].z()*UMean_[n].z());

                pField_[n] = p[n];

                scalarPressure_[n] = (1.0/3.0)*(p[n].xx() + p[n].yy() + p[n].zz());

                // make reference
                tensorField tau(uspMols.size(), Zero);

                tau[n] = -p[n];
                tau[n].xx() += scalarPressure_[n];
                tau[n].yy() += scalarPressure_[n];
                tau[n].zz() += scalarPressure_[n];
                tauField_[n] = tau[n];

                vectorField q(mols_.size(), Zero);

                q[n].x() =
                    rhoN_[n]
                   *(
                        0.5*(mccu[n]/mols[n])
                      - 0.5*(mcc[n]/mols[n])*UMean_[n].x()
                      + eu[n]/mols[n]
                      - (e[n]/mols[n])*UMean_[n].x()
                    )
                  - p[n].xx()*UMean_[n].x()
                  - p[n].xy()*UMean_[n].y()
                  - p[n].xz()*UMean_[n].z();

                 //terms involving pressure tensor should not be multiplied by the number density (see Bird corrigendum)


                q[n].y() =
                    rhoN_[n]
                  *(
                       0.5*(mccv[n]/(mols[n]))
                     - 0.5*(mcc[n]/(mols[n]))*UMean_[n].y()
                     + ev[n]/(mols[n])
                     - (e[n]/(mols[n]))*UMean_[n].y()
                    )
                  - p[n].yx()*UMean_[n].x()
                  - p[n].yy()*UMean_[n].y()
                  - p[n].yz()*UMean_[n].z();

                q[n].z() =
                    rhoN_[n]
                   *(
                        0.5*(mccw[n]/(mols[n]))
                      - 0.5*(mcc[n]/(mols[n]))*UMean_[n].z()
                      + ew[n]/(mols[n])
                      - (e[n]/(mols[n]))*UMean_[n].z()
                    )
                  - p[n].zx()*UMean_[n].x()
                  - p[n].zy()*UMean_[n].y()
                  - p[n].zz()*UMean_[n].z();

                qField_[n] = q[n];

                vectorField qInternal(mols_.size(), Zero);

                qInternal[n].x() =
                    rhoN_[n]*(eu[n]/(mols[n]) - (e[n]/(mols[n]))*UMean_[n].x());

                qInternal[n].y() =
                    rhoN_[n]*(ev[n]/(mols[n]) - (e[n]/(mols[n]))*UMean_[n].y());

                qInternal[n].z() =
                    rhoN_[n]*(ew[n]/(mols[n]) - (e[n]/(mols[n]))*UMean_[n].z());

                qInternalField_[n] = qInternal[n];


                vectorField qTranslational(mols_.size(), Zero);

                qTranslational[n].x() =
                    rhoN_[n]
                   *(
                        0.5*(mccu[n]/(mols[n]))
                      - 0.5*(mcc[n]/(mols[n]))*UMean_[n].x()
                    )
                  - p[n].xx()*UMean_[n].x()
                  - p[n].xy()*UMean_[n].y()
                  - p[n].xz()*UMean_[n].z();

                qTranslational[n].y() =
                    rhoN_[n]
                   *(
                        0.5*(mccv[n]/(mols[n]))
                      - 0.5*(mcc[n]/(mols[n]))*UMean_[n].y()
                    )
                  - p[n].yx()*UMean_[n].x()
                  - p[n].yy()*UMean_[n].y()
                  - p[n].yz()*UMean_[n].z();

                qTranslational[n].z() =
                    rhoN_[n]
                   *(
                        0.5*(mccw[n]/(mols[n]))
                      - 0.5*(mcc[n]/(mols[n]))*UMean_[n].z()
                    )
                  - p[n].zx()*UMean_[n].x()
                  - p[n].zy()*UMean_[n].y()
                  - p[n].zz()*UMean_[n].z();

                qTranslationalField_[n] = qTranslational[n];

                const label nBins = binModel_->nBins();

                // vibrational temperature
                scalarList totalvDof(nBins, Zero);
                scalarList totalvDofOverall(nBins, Zero);
                scalarList vibT(nBins, Zero);
                List<scalarList> degreesOfFreedomSpecies;
                List<scalarList> vibTID;

                degreesOfFreedomSpecies.setSize(nBins);
                vibTID.setSize(nBins);

                forAll(degreesOfFreedomSpecies, b)
                {
                    degreesOfFreedomSpecies[b].setSize(typeIds_.size(), Zero);
                    vibTID[b].setSize(typeIds_.size(), Zero);
                }

                List<List<scalarList>> degreesOfFreedomMode;
                List<List<scalarList>> vibTMode;

                degreesOfFreedomMode.setSize(nBins);
                vibTMode.setSize(nBins);

                forAll(degreesOfFreedomMode, b)
                {
                    degreesOfFreedomMode[b].setSize(typeIds_.size());
                    vibTMode[b].setSize(typeIds_.size());

                    forAll(degreesOfFreedomMode[b], iD)
                    {
                        degreesOfFreedomMode[b][iD].setSize
                        (
                            vibrationalETotal.size(), Zero
                        );
                        vibTMode[b][iD].setSize(vibrationalETotal.size(), Zero);
                    }
                }

                forAll(vibrationalETotal[n], iD)
                {
                    forAll(vibrationalETotal[n][iD], m)
                    {
                        if
                        (
                            vibrationalETotal[n][iD][m] > VSMALL
                         && speciesMols[n][iD] > VSMALL
                        )
                        {
                            scalar thetaV =
                                cloud_.constProps(typeIds_[iD]).thetaV()[m];

                            scalarList vibrationalEMean =
                                (vibrationalETotal[n][iD]/speciesMols[n][iD]);

                            scalar iMean = 0.0;

                            iMean =
                                vibrationalEMean[m]
                               /(physicoChemical::k.value()*thetaV);

                            vibTMode[n][iD][m] =
                                thetaV / log(1.0 + (1.0/iMean));

                            degreesOfFreedomMode[n][iD][m] =
                                (2.0*thetaV/vibTMode[n][iD][m])
                               /(exp(thetaV/vibTMode[n][iD][m]) - 1.0);
                        }
                    }

                    forAll(vibrationalETotal[n][iD], m)
                    {
                        degreesOfFreedomSpecies[n][iD] +=
                            degreesOfFreedomMode[n][iD][m];
                    }

                    forAll(vibrationalETotal[n][iD], m)
                    {
                        if (degreesOfFreedomSpecies[n][iD] > VSMALL)
                        {
                            vibTID[n][iD] +=
                                vibTMode[n][iD][m]
                               *degreesOfFreedomMode[n][iD][m]
                               /degreesOfFreedomSpecies[n][iD];
                        }
                    }

                    totalvDof[n] += degreesOfFreedomSpecies[n][iD];

                    scalar fraction = 0.0;

                    if (molsInt[n] > VSMALL)
                    {
                        fraction = speciesMols[n][iD]/molsInt[n];
                    }

                    scalar fractionOverall = speciesMols[n][iD]/mols[n];

                    if (fraction > SMALL)
                    {
                        totalvDofOverall[n] +=
                            totalvDof[n]*(fractionOverall/fraction);
                    }

                    vibT[n] += vibTID[n][iD]*fraction;
                }

                vibrationalTemperature_[n] = vibT[n];

                // electronic temperature
                scalarList totalEDof(nBins, Zero);
                scalarList elecT(nBins, Zero);

                forAll(speciesMols[n], iD)
                {
                    label nElectronicLevels =
                        cloud_.constProps(typeIds_[iD]).nElectronicLevels();

                    if
                    (
                        nElectronicLevels > 1
                     && speciesMols[n][iD] > VSMALL
                     && molsElec[n] > VSMALL
                    )
                    {
                        const scalarList& electronicEnergies =
                            cloud_.constProps(typeIds_[iD]).
                                electronicEnergyList();
                        const labelList& degeneracies =
                            cloud_.constProps(typeIds_[iD]).degeneracyList();

                        if
                        (
                            nParticlesGroundElectronicState[n][iD] > VSMALL
                         && nParticlesFirstElectronicState[n][iD] > VSMALL
                         && (
                                (nParticlesGroundElectronicState[n][iD]*degeneracies[1]) != (nParticlesFirstElectronicState[n][iD]*degeneracies[0]))
                            )
                        {
                            scalar fraction = speciesMols[n][iD]/molsElec[n];

                            scalar elecTID =
                                (electronicEnergies[1]-electronicEnergies[0])
                               /(
                                   physicoChemical::k.value()
                                  *log
                                   (
                                       (nParticlesGroundElectronicState[n][iD]*degeneracies[1])
                                      /(nParticlesFirstElectronicState[n][iD]*degeneracies[0])
                                    )
                                );

                            if (elecTID > VSMALL)
                            {
                                elecT[n] += fraction*elecTID;
                            }

                            scalar eDof =
                                (2.0*(electronicETotal[n][iD]/speciesMols[n][iD]))
                               /(physicoChemical::k.value()*elecTID);

                            totalEDof[n] += fraction*eDof;
                        }
                    }
                }

                electronicTemperature_[n] = elecT[n];


                //overallTemperature

                scalar nRotDof = 0.0;

                if (uspMols[n] > VSMALL)
                {
                    nRotDof = rotationalDofMean[n]/uspMols[n];
                }

                overallTemperature_[n] =
                    (
                        3.0*translationalTemperature_[n]
                      + nRotDof*rotationalTemperature_[n]
                      + totalvDof[n]*vibrationalTemperature_[n]
                      + totalEDof[n]*electronicTemperature_[n]
                    )
                   /(3.0 + nRotDof + totalvDof[n] + totalEDof[n]);

                forAll(mfp_[n], iD)
                {
                    const auto& cpi = cloud_.constProps(typeIds_[iD]);

                    for (label qspec=0; qspec<typeIds_.size(); ++qspec)
                    {
                        const auto& cpq = cloud_.constProps(qspec);

                        scalar dPQ = 0.5*(cpi.d() + cpq.d());
                        scalar omegaPQ = 0.5*(cpi.omega() + cpq.omega());
                        scalar massRatio = cpi.mass()/cpq.mass();

                        if
                        (
                            speciesMols[n][qspec] > VSMALL
                         && translationalTemperature_[n] > VSMALL
                        )
                        {
                            scalar nDensQ =
                                cloud_.nParticle()*speciesMols[n][qspec]
                               /(volume*averagingCounter_);

                            //Bird, eq (4.76)
                            mfp_[n][iD] +=
                                pi*dPQ*dPQ*nDensQ
                               *pow
                               (
                                   273.0/translationalTemperature_[n],
                                   omegaPQ - 0.5
                                )
                               *sqrt(1.0 + massRatio);
                        }
                    }

                    if (mfp_[n][iD] > VSMALL)
                    {
                        mfp_[n][iD] = 1.0/mfp_[n][iD];
                    }
                }

                meanFreePath_[n] = 0.0;

                forAll(mfp_[n], iD)
                {
                    if(rhoN_[n] > VSMALL)
                    {
                        scalar nDensP =
                            cloud_.nParticle()*speciesMols[n][iD]
                           /(volume*averagingCounter_);

                        // Bird, eq (4.77)
                        meanFreePath_[n] += mfp_[n][iD]*nDensP/rhoN_[n];
                    }
                    else
                    {
                        meanFreePath_[n] = GREAT;
                    }
                }

                mfp_[n] = scalar(0.0);

                scalar molecularMass = 0.0;
                scalar molarconstantPressureSpecificHeat = 0.0;
                scalar molarconstantVolumeSpecificHeat = 0.0;
                scalar gasConstant = 0.0;
                scalar gamma = 0.0;
                scalar speedOfSound = 0.0;

                forAll(mfp_[n], iD)
                {
                    if (rhoN_[n] > VSMALL)
                    {
                        const auto& cpi = cloud_.constProps(typeIds_[iD]);

                        molecularMass +=
                            cpi.mass()*speciesMols[n][iD]/uspMols[n];
                        molarconstantPressureSpecificHeat +=
                            (5.0 + cpi.rotationalDoF())
                           *(speciesMols[n][iD]/uspMols[n]);
                        molarconstantVolumeSpecificHeat +=
                            (3.0 + cpi.rotationalDoF())
                           *(speciesMols[n][iD]/uspMols[n]);
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
        }

        if (timeVel_.resetFieldsAtOutput() && (mesh_.time().value() < timeVel_.resetFieldsAtOutputUntilTime()+0.5*mesh_.time().deltaTValue()))
        {
            // Reset fields
            averagingCounter_ = 0.0;

            volumeOfCellsInBin_ = 0.0;
            mols_ = 0.0;
            uspMols_ = 0.0;
            molsInt_ = 0.0;
            molsElec_ = 0.0;
            mass_ = 0.0;
            mcc_ = 0.0;
            mom_ = vector::zero;
            UCollected_ = vector::zero;
            rotationalEMean_ = 0.0;
            rotationalDofMean_ = 0.0;

            muu_ = 0.0;
            muv_ = 0.0;
            muw_ = 0.0;
            mvv_ = 0.0;
            mvw_ = 0.0;
            mww_ = 0.0;
            mccu_ = 0.0;
            mccv_ = 0.0;
            mccw_ = 0.0;
            mcc_ = 0.0;
            eu_ = 0.0;
            ev_ = 0.0;
            ew_ = 0.0;
            e_ = 0.0;

            forAll(electronicETotal_, b)
            {
                electronicETotal_[b] = 0.0;
                nParticlesGroundElectronicState_[b] = 0.0;
                nParticlesFirstElectronicState_[b] = 0.0;
                speciesMols_[b] = 0.0;
                mccSpecies_[b] = 0.0;

                forAll(vibrationalETotal_[b], iD)
                {
                    vibrationalETotal_[b][iD] = 0.0;
                }
            }
        }

    }

    if (averagingAcrossManyRuns_)
    {
        writeOut();
    }

}


void Foam::uspBinsMethod::writeField()
{
    const Time& runTime = mesh_.time();

    if (runTime.writeTime() && Pstream::master())
    {
        scalarField bins(binModel_->binPositions());
        vectorField vectorBins(binModel_->bins());

        const word prefix("bins_OneDim_" + regionName_ + "_" + fieldName_);

        // Output densities
        if (outputField_[0])
        {
            writeTimeData
            (
                timePath_,
                prefix+"_N.xy",
                bins,
                N_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_N_3D_pos.xy",
                vectorBins,
                N_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_rhoN.xy",
                bins,
                rhoN_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_rhoN_3D_pos.xy",
                vectorBins,
                rhoN_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_rhoM.xy",
                bins,
                rhoM_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_rhoM_3D_pos.xy",
                vectorBins,
                rhoM_
            );
        }

        // Output velocities
        if (outputField_[1])
        {
            writeTimeData
            (
                timePath_,
                prefix+"_U_SAM.xyz",
                bins,
                UMean_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_U_SAM_3D_pos.xyz",
                vectorBins,
                UMean_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_U_CAM.xyz",
                bins,
                UCAM_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_U_CAM_3D_pos.xyz",
                vectorBins,
                UCAM_
            );
        }

        // Output temperature
        if (outputField_[2])
        {
            writeTimeData
            (
                timePath_,
                prefix+"_translationalTemperature.xy",
                bins,
                translationalTemperature_
            );
            writeTimeData
            (
                timePath_,
                prefix+"_translationalTemperature_3D_pos.xyz",
                vectorBins,
                translationalTemperature_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_rotationalTemperature.xy",
                bins,
                rotationalTemperature_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_rotationalTemperature_3D_pos.xyz",
                vectorBins,
                rotationalTemperature_
            );
            writeTimeData
            (
                timePath_,
                prefix+"_vibrationalTemperature.xy",
                bins,
                vibrationalTemperature_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_vibrationalTemperature_3D_pos.xyz",
                vectorBins,
                vibrationalTemperature_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_electronicTemperature.xy",
                bins,
                electronicTemperature_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_electronicTemperature_3D_pos.xyz",
                vectorBins,
                electronicTemperature_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_overallTemperature.xy",
                bins,
                overallTemperature_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_overallTemperature_3D_pos.xyz",
                vectorBins,
                overallTemperature_
            );
        }

        // output pressure
        if (outputField_[3])
        {
            writeTimeData
            (
                timePath_,
                prefix+"_pressureTensor.xyz",
                bins,
                pField_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_pressureTensor_3D_pos.xyz",
                vectorBins,
                pField_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_p.xy",
                bins,
                scalarPressure_
            );
            writeTimeData
            (
                timePath_,
                prefix+"_p_3D_pos.xyz",
                vectorBins,
                scalarPressure_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_stressTensor.xyz",
                bins,
                tauField_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_stressTensor_3D_pos.xyz",
                vectorBins,
                tauField_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_heatFluxVector.xyz",
                bins,
                qField_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_heatFluxVector_3D_pos.xyz",
                vectorBins,
                qField_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_internalHeatFluxVector.xyz",
                bins,
                qInternalField_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_internalHeatFluxVector_3D_pos.xyz",
                vectorBins,
                qInternalField_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_translationalHeatFluxVector.xyz",
                bins,
                qTranslationalField_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_translationalHeatFluxVector_3D_pos.xyz",
                vectorBins,
                qTranslationalField_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_variableHardSphereMeanFreePath.xy",
                bins,
                meanFreePath_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_variableHardSphereMeanFreePath_3D_pos.xyz",
                vectorBins,
                meanFreePath_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_Ma.xy",
                bins,
                Ma_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_Ma_3D_pos.xyz",
                vectorBins,
                Ma_
            );
        }

        if (permeabilityMeasurements_)
        {
            writeTimeData
            (
                timePath_,
                prefix+"_binVolume.xy",
                bins,
                binVolume_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_binVolume_3D_pos.xyz",
                vectorBins,
                binVolume_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_binNumberOfParticles.xy",
                bins,
                N_*cloud_.nParticle()
            );

            writeTimeData
            (
                timePath_,
                prefix+"_binNumberOfParticles_3D_pos.xyz",
                vectorBins,
                N_*cloud_.nParticle()
            );

            writeTimeData
            (
                timePath_,
                prefix+"_binMass.xy",
                bins,
                averageMass_
            );

            writeTimeData
            (
                timePath_,
                prefix+"_binMass_3D_pos.xyz",
                vectorBins,
                averageMass_
            );
        }
    }
}


void Foam::uspBinsMethod::updateProperties(const dictionary& dict)
{
    // The main properties should be updated first
    uspField::updateProperties(dict);
}


// ************************************************************************* //
