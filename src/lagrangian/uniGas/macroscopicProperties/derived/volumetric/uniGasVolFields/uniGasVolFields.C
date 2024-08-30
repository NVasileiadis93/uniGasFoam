/*---------------------------------------------------------------------------* \
  =========                 |
  \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration     |
    \  /    A nd           | www.openfoam.com
     \/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "uniGasVolFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasVolFields, 0);

addToRunTimeSelectionTable(uniGasField, uniGasVolFields, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasVolFields::uniGasVolFields
(
    const Time& t,
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.get<word>("field")),
    typeIds_(cloud_.getTypeIDs(propsDict_)),
    Tref_(propsDict_.get<scalar>("Tref")),
    densityOnly_(propsDict_.getOrDefault<bool>("densityOnly",false)),
    measureMeanFreePath_(propsDict_.getOrDefault<bool>("measureMeanFreePath",false)),
    measureErrors_(propsDict_.getOrDefault<bool>("measureErrors",false)),
    averagingAcrossManyRuns_(propsDict_.getOrDefault<bool>("averagingAcrossManyRuns",false)),
    sampleCounter_(0),
    nAvTimeSteps_(0),
    timeAvCounter_(0.0),
    n_(),
    t1_(),
    t2_(),
    uniGasRhoNMean_
    (
        IOobject
        (
            "uniGasRhoNMean_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimVolume, Zero)
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN_" + fieldName_,
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
            "rhoM_" + fieldName_,
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
            "p_" + fieldName_,
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
            "translationalT_" + fieldName_,
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
            "rotationalT_" + fieldName_,
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
            "vibrationalT_" + fieldName_,
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
            "electronicT_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, Zero)
    ),
    overallT_
    (
        IOobject
        (
            "overallT_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTemperature, Zero)
    ),
    q_
    (
        IOobject
        (
            "surfaceHeatTransfer_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), Zero)
    ),
    tau_
    (
        IOobject
        (
            "surfaceShearStress_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure, Zero)
    ),
    MFP_
    (
        IOobject
        (
            "variableHardSphereMeanFreePath_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimLength, Zero)
    ),
    dxMFP_
    (
        IOobject
        (
            "subCellSizeMFPRatio_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    MCR_
    (
        IOobject
        (
            "meanCollisionRate_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    ),
    MCT_
    (
        IOobject
        (
            "meanCollisionTime_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimTime, Zero)
    ),
    dtMCT_
    (
        IOobject
        (
            "timeStepMCTRatio_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    Ma_
    (
        IOobject
        (
            "Ma_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    densityError_
    (
        IOobject
        (
            "densityError_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    velocityError_
    (
        IOobject
        (
            "velocityError_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    temperatureError_
    (
        IOobject
        (
            "temperatureError_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    pressureError_
    (
        IOobject
        (
            "pressureError_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    UMean_
    (
        IOobject
        (
            "UMean_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimLength/dimTime, Zero)
    ),
    fD_
    (
        IOobject
        (
            "fD_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimensionSet(1, -1, -2, 0, 0), Zero)
    ),
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
    momentumMean_(mesh.nCells(), Zero),
    momentumMeanXnParticle_(mesh.nCells(), Zero),
    boundaryCells_(),
    vibrationalETotal_(),
    electronicETotal_(),
    nParcels_(),
    nParcelsXnParticle_(),
    mccSpecies_(),
    nGroundElectronicLevel_(),
    nFirstElectronicLevel_(),
    mfp_(),
    mcr_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    rotationalEBF_(),
    rotationalDofBF_(),
    qBF_(),
    totalvDofBF_(),
    speciesRhoNIntBF_(),
    speciesRhoNElecBF_(),
    momentumBF_(),
    fDBF_(),
    vibrationalEBF_(),
    electronicEBF_(),
    speciesRhoNBF_(),
    mccSpeciesBF_()
{

    // Note; outer list is typeIds, inner list is number of cells on the
    // mesh
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

    mfp_.setSize(typeIds_.size());

    for (auto& m : mfp_)
    {
        m.setSize(mesh_.nCells());
    }

    mcr_.setSize(typeIds_.size());

    for (auto& m : mcr_)
    {
        m.setSize(mesh_.nCells());
    }

    boundaryCells_.setSize(mesh.boundaryMesh().size());

    forAll(boundaryCells_, p)
    {
        const polyPatch& patch = mesh.boundaryMesh()[p];

        boundaryCells_[p].setSize(patch.size());

        forAll(boundaryCells_[p], c)
        {
            boundaryCells_[p][c] = patch.faceCells()[c];
        }
    }

    // initialisation
    rhoNBF_.setSize(mesh_.boundaryMesh().size());
    rhoMBF_.setSize(mesh_.boundaryMesh().size());
    linearKEBF_.setSize(mesh_.boundaryMesh().size());
    momentumBF_.setSize(mesh_.boundaryMesh().size());
    rotationalEBF_.setSize(mesh_.boundaryMesh().size());
    rotationalDofBF_.setSize(mesh_.boundaryMesh().size());
    qBF_.setSize(mesh_.boundaryMesh().size());
    fDBF_.setSize(mesh_.boundaryMesh().size());
    totalvDofBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNIntBF_.setSize(mesh_.boundaryMesh().size());
    speciesRhoNElecBF_.setSize(mesh_.boundaryMesh().size());
    n_.setSize(mesh_.boundaryMesh().size());
    t1_.setSize(mesh_.boundaryMesh().size());
    t2_.setSize(mesh_.boundaryMesh().size());

    forAll(rhoNBF_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];

        rhoNBF_[j].setSize(patch.size(), 0.0);
        rhoMBF_[j].setSize(patch.size(), 0.0);
        linearKEBF_[j].setSize(patch.size(), 0.0);
        momentumBF_[j].setSize(patch.size(), Zero);
        rotationalEBF_[j].setSize(patch.size(), 0.0);
        rotationalDofBF_[j].setSize(patch.size(), 0.0);
        qBF_[j].setSize(patch.size(), 0.0);
        fDBF_[j].setSize(patch.size(), Zero);
        totalvDofBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNIntBF_[j].setSize(patch.size(), 0.0);
        speciesRhoNElecBF_[j].setSize(patch.size(), 0.0);
        n_[j].setSize(patch.size(), Zero);
        t1_[j].setSize(patch.size(), Zero);
        t2_[j].setSize(patch.size(), Zero);
    }

    calculateWallUnitVectors();

    vibrationalEBF_.setSize(typeIds_.size());
    electronicEBF_.setSize(typeIds_.size());
    speciesRhoNBF_.setSize(typeIds_.size());
    mccSpeciesBF_.setSize(typeIds_.size());

    forAll(vibrationalEBF_, i)
    {
        vibrationalEBF_[i].setSize(mesh_.boundaryMesh().size());
        electronicEBF_[i].setSize(mesh_.boundaryMesh().size());
        speciesRhoNBF_[i].setSize(mesh_.boundaryMesh().size());
        mccSpeciesBF_[i].setSize(mesh_.boundaryMesh().size());

        forAll(vibrationalEBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            vibrationalEBF_[i][j].setSize(patch.size(), 0.0);
            electronicEBF_[i][j].setSize(patch.size(), 0.0);
            speciesRhoNBF_[i][j].setSize(patch.size(), 0.0);
            mccSpeciesBF_[i][j].setSize(patch.size(), 0.0);
        }
    }

    // read in stored data from dictionary
    if (averagingAcrossManyRuns_)
    {
        readIn();
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasVolFields::readIn()
{

    localIOdictionary dict
    (
        IOobject
        (
            "volFieldsMethod_" +fieldName_,
            mesh_.time().timeName(),
            "uniform",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    dict.readIfPresent("nTimeSteps", nAvTimeSteps_);
    dict.readIfPresent("timeCounter", timeAvCounter_);
    dict.readIfPresent("rhoNMean", rhoNMean_);
    dict.readIfPresent("rhoNMeanXnParticle", rhoNMeanXnParticle_);
    dict.readIfPresent("rhoNMeanInt", rhoNMeanInt_);
    dict.readIfPresent("molsElec", molsElec_);
    dict.readIfPresent("rhoMMean", rhoMMean_);
    dict.readIfPresent("rhoMMeanXnParticle", rhoMMeanXnParticle_);
    dict.readIfPresent("linearKEMean", linearKEMean_);
    dict.readIfPresent("linearKEMeanXnParticle", linearKEMeanXnParticle_);
    dict.readIfPresent("rotationalEMean", rotationalEMean_);
    dict.readIfPresent("rotationalDofMean", rotationalDofMean_);
    dict.readIfPresent("momentumMean", momentumMean_);
    dict.readIfPresent("momentumMeanXnParticle", momentumMeanXnParticle_);
    dict.readIfPresent("vibrationalETotal", vibrationalETotal_);
    dict.readIfPresent("electronicETotal", electronicETotal_);
    dict.readIfPresent("nParcels", nParcels_);
    dict.readIfPresent("nParcelsXnParticle", nParcelsXnParticle_);
    dict.readIfPresent("mccSpecies", mccSpecies_);
    dict.readIfPresent("nGroundElectronicLevel", nGroundElectronicLevel_);
    dict.readIfPresent("nFirstElectronicLevel", nFirstElectronicLevel_);
    dict.readIfPresent("mfp", mfp_);
    dict.readIfPresent("mcr", mcr_);
    dict.readIfPresent("rhoNBF", rhoNBF_);
    dict.readIfPresent("rhoMBF", rhoMBF_);
    dict.readIfPresent("linearKEBF", linearKEBF_);
    dict.readIfPresent("rotationalEBF", rotationalEBF_);
    dict.readIfPresent("rotationalDofBF", rotationalDofBF_);
    dict.readIfPresent("qBF", qBF_);
    dict.readIfPresent("totalvDofBF", totalvDofBF_);
    dict.readIfPresent("speciesRhoNIntBF", speciesRhoNIntBF_);
    dict.readIfPresent("speciesRhoNElecBF", speciesRhoNElecBF_);
    dict.readIfPresent("momentumBF", momentumBF_);
    dict.readIfPresent("fDBF", fDBF_);
    dict.readIfPresent("vibrationalEBF", vibrationalEBF_);
    dict.readIfPresent("electronicEBF", electronicEBF_);
    dict.readIfPresent("speciesRhoNBF", speciesRhoNBF_);
    dict.readIfPresent("mccSpeciesBF", mccSpeciesBF_);
}


void Foam::uniGasVolFields::writeOut()
{
    if (mesh_.time().writeTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "volFieldsMethod_" +fieldName_,
                mesh_.time().timeName(),
                "uniform",
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );

        dict.add("nTimeSteps", nAvTimeSteps_);
        dict.add("timeCounter", timeAvCounter_);
        dict.add("rhoNMean", rhoNMean_);
        dict.add("rhoNMeanXnParticle", rhoNMeanXnParticle_);
        dict.add("rhoNMeanInt", rhoNMeanInt_);
        dict.add("molsElec", molsElec_);
        dict.add("rhoMMean", rhoMMean_);
        dict.add("rhoMMeanXnParticle", rhoMMeanXnParticle_);
        dict.add("linearKEMean", linearKEMean_);
        dict.add("linearKEMeanXnParticle", linearKEMeanXnParticle_);
        dict.add("rotationalEMean", rotationalEMean_);
        dict.add("rotationalDofMean", rotationalDofMean_);
        dict.add("momentumMean", momentumMean_);
        dict.add("momentumMeanXnParticle", momentumMeanXnParticle_);
        dict.add("vibrationalETotal", vibrationalETotal_);
        dict.add("electronicETotal", electronicETotal_);
        dict.add("nParcels", nParcels_);
        dict.add("nParcelsXnParticle", nParcelsXnParticle_);
        dict.add("mccSpecies", mccSpecies_);
        dict.add("nGroundElectronicLevel", nGroundElectronicLevel_);
        dict.add("nFirstElectronicLevel", nFirstElectronicLevel_);
        dict.add("mfp", mfp_);
        dict.add("mcr", mcr_);
        dict.add("rhoNBF", rhoNBF_);
        dict.add("rhoMBF", rhoMBF_);
        dict.add("linearKEBF", linearKEBF_);
        dict.add("rotationalEBF", rotationalEBF_);
        dict.add("rotationalDofBF", rotationalDofBF_);
        dict.add("qBF", qBF_);
        dict.add("totalvDofBF", totalvDofBF_);
        dict.add("speciesRhoNIntBF", speciesRhoNIntBF_);
        dict.add("speciesRhoNElecBF", speciesRhoNElecBF_);
        dict.add("momentumBF", momentumBF_);
        dict.add("fDBF", fDBF_);
        dict.add("vibrationalEBF", vibrationalEBF_);
        dict.add("electronicEBF", electronicEBF_);
        dict.add("speciesRhoNBF", speciesRhoNBF_);
        dict.add("mccSpeciesBF", mccSpeciesBF_);

        dict.regIOobject::writeObject
        (
            IOstreamOption(mesh_.time().writeFormat()),
            true
        );
    }
}


void Foam::uniGasVolFields::calculateWallUnitVectors()
{
    forAll(n_, patchi)
    {
        const polyPatch& pPatch = mesh_.boundaryMesh()[patchi];

        if (isA<wallPolyPatch>(pPatch))
        {
            const vectorField& fC = pPatch.faceCentres();

            forAll(n_[patchi], facei)
            {
                n_[patchi][facei] = pPatch.faceAreas()[facei];
                n_[patchi][facei].normalise();

                // Wall tangential unit vector. Use the direction between the
                // face centre and the first vertex in the list
                t1_[patchi][facei] = fC[facei] -
                    mesh_.points()[mesh_.faces()[pPatch.start() + facei][0]];
                t1_[patchi][facei].normalise();

                //  Other tangential unit vector.  Rescaling in case face is not
                //  flat and n and t1 aren't perfectly orthogonal
                t2_[patchi][facei] = n_[patchi][facei]^t1_[patchi][facei];
                t2_[patchi][facei].normalise();
            }
        }
    }
}


void Foam::uniGasVolFields::createField()
{
    
    Info << "Initialising uniGasVolFields field" << endl;

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


void Foam::uniGasVolFields::calculateField()
{
    sampleCounter_++;

    const scalar& deltaT = mesh_.time().deltaTValue();

    if (sampleInterval_ <= sampleCounter_)
    {

        nAvTimeSteps_++;

        timeAvCounter_ += deltaT;

        //obtain cell measurements
        auto& cm = cloud_.cellPropMeasurements();

        if (densityOnly_)
        {

            forAll(cm.rhoNMean(), i)
            {

                const label iD = typeIds_.find(i);
                
                if (iD != -1)
                {

                    rhoNMean_ += deltaT*cm.rhoNMean()[i];
                    rhoNMeanXnParticle_ += deltaT*cm.rhoNMeanXnParticle()[i];
                    rhoMMeanXnParticle_ += deltaT*cm.rhoMMeanXnParticle()[i];
                }

            }

        }
        else
        {

            forAll(cm.rhoNMean(), i)
            {

                const label iD = typeIds_.find(i);

                if (iD != -1)
                {

                    rhoNMean_ += deltaT*cm.rhoNMean()[i];
                    rhoMMean_ += deltaT*cm.rhoMMean()[i];
                    linearKEMean_ += deltaT*cm.linearKEMean()[i];
                    momentumMean_ += deltaT*cm.momentumMean()[i];
                    rotationalEMean_ += deltaT*cm.rotationalEMean()[i];
                    rotationalDofMean_ += deltaT*cm.rotationalDofMean()[i];
                    electronicETotal_[iD] += deltaT*cm.electronicETotal()[i];
                    rhoNMeanXnParticle_ += deltaT*cm.rhoNMeanXnParticle()[i];
                    rhoMMeanXnParticle_ += deltaT*cm.rhoMMeanXnParticle()[i];
                    momentumMeanXnParticle_ += deltaT*cm.momentumMeanXnParticle()[i];
                    linearKEMeanXnParticle_ += deltaT*cm.linearKEMeanXnParticle()[i];

                    rhoNMeanInt_ += deltaT*cm.rhoNMeanInt()[i];
                    molsElec_ += deltaT*cm.molsElec()[i];

                    nParcels_[iD] += deltaT*cm.nParcels()[i];
                    nParcelsXnParticle_[iD] += deltaT*cm.nParcelsXnParticle()[i];
                    mccSpecies_[iD] += deltaT*cm.mccSpecies()[i];

                    nGroundElectronicLevel_[iD] += deltaT*cm.nGroundElectronicLevel()[i];
                    nFirstElectronicLevel_[iD] += deltaT*cm.nFirstElectronicLevel()[i];

                    forAll(vibrationalETotal_[iD], v)
                    {
                        vibrationalETotal_[iD][v] += deltaT*cm.vibrationalETotal()[i][v];
                    }
                }

            }

        }

        // obtain boundary measurements
        auto& bm = cloud_.boundaryFluxMeasurements();

        forAll(bm.rhoNBF(), i)
        {

            const label iD = typeIds_.find(i);

            if (iD != -1)
            {

                forAll(bm.rhoNBF()[i], j)
                {
                    forAll(bm.rhoNBF()[i][j], k)
                    {
                        rhoNBF_[j][k] += deltaT*bm.rhoNBF()[i][j][k];
                        rhoMBF_[j][k] += deltaT*bm.rhoMBF()[i][j][k];
                        linearKEBF_[j][k] += deltaT*bm.linearKEBF()[i][j][k];
                        momentumBF_[j][k] += deltaT*bm.momentumBF()[i][j][k];
                        rotationalEBF_[j][k] += deltaT*bm.rotationalEBF()[i][j][k];
                        rotationalDofBF_[j][k] += deltaT*bm.rotationalDofBF()[i][j][k];
                        qBF_[j][k] += deltaT*bm.qBF()[i][j][k];
                        fDBF_[j][k] += deltaT*bm.fDBF()[i][j][k];
                        speciesRhoNBF_[iD][j][k] += deltaT*bm.rhoNBF()[i][j][k];
                        vibrationalEBF_[iD][j][k] += deltaT*bm.vibrationalEBF()[i][j][k];
                        electronicEBF_[iD][j][k] += deltaT*bm.electronicEBF()[i][j][k];
                        mccSpeciesBF_[iD][j][k] += deltaT*bm.mccSpeciesBF()[i][j][k];
                        speciesRhoNIntBF_[j][k] += deltaT*bm.rhoNIntBF()[i][j][k];
                        speciesRhoNElecBF_[j][k] += deltaT*bm.rhoNElecBF()[i][j][k];
                    }
                }

            }
        }

        sampleCounter_ = 0;
    }

    if (mesh_.time().writeTime())
    {

        if (densityOnly_)
        {
            forAll(rhoNMean_, cell)
            {
                if (rhoNMean_[cell] > VSMALL)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[cell];

                    uniGasRhoNMean_[cell] = rhoNMean_[cell]/timeAvCounter_;

                    rhoN_[cell] = (rhoNMeanXnParticle_[cell])/(timeAvCounter_*cellVolume);

                    rhoM_[cell] = (rhoMMeanXnParticle_[cell])/(timeAvCounter_*cellVolume);

                }
                else
                {
                    uniGasRhoNMean_[cell] = 0.0;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                }

            }
        }
        else
        {

            forAll(rhoNMean_, cell)
            {
                if (rhoNMean_[cell] > VSMALL)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[cell];

                    uniGasRhoNMean_[cell] = rhoNMean_[cell]/timeAvCounter_;

                    rhoN_[cell] =
                        (rhoNMeanXnParticle_[cell])/(timeAvCounter_*cellVolume);

                    rhoM_[cell] =
                        (rhoMMeanXnParticle_[cell])/(timeAvCounter_*cellVolume);

                    scalar rhoMMean =
                        rhoMMeanXnParticle_[cell]/(cellVolume*timeAvCounter_);
                    UMean_[cell] =
                        momentumMeanXnParticle_[cell]
                       /(rhoMMean*cellVolume*timeAvCounter_);
                    scalar linearKEMean =
                        0.5*linearKEMeanXnParticle_[cell]
                       /(cellVolume*timeAvCounter_);
                    scalar rhoNMean =
                        rhoNMeanXnParticle_[cell]/(cellVolume*timeAvCounter_);

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
                    // not zero so that weighted decomposition still works
                    uniGasRhoNMean_[cell] = 0.001;
                    rhoN_[cell] = 0.0;
                    rhoM_[cell] = 0.0;
                    UMean_[cell] = vector::zero;
                    translationalT_[cell] = 0.0;
                    p_[cell] = 0.0;
                }

                // Rotational temperature
                if (rotationalDofMean_[cell] > VSMALL && timeAvCounter_ > VSMALL)
                {
                    scalar rotationalEMean = rotationalEMean_[cell]/timeAvCounter_;
                    
                    scalar rotationalDofMean = rotationalDofMean_[cell]/timeAvCounter_;

                    rotationalT_[cell] = (2.0/physicoChemical::k.value())*(rotationalEMean/rotationalDofMean);
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
                        vibT += vibTID[iD]*nParcels_[iD][cell]/rhoNMeanInt_[cell];
                    }
                }

                vibrationalT_[cell] = vibT;             

                // Electronic temperature
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

                // Overall temperature
                overallT_[cell] =
                    (
                        (3.0*translationalT_[cell])
                      + (nRotDof*rotationalT_[cell])
                      + (totalvDof*vibrationalT_[cell])
                      + (totalEDof*electronicT_[cell])
                    )
                   /(3.0 + nRotDof + totalvDof + totalEDof);

                // Mach number
                scalar Cp = 0.0;
                scalar Cv = 0.0;
                scalar molecularMass = 0.0;
                scalar Cv_p = 0.0;

                forAll(nParcels_, iD)
                {

                    if (rhoNMean_[cell] > VSMALL)
                    {
                        molecularMass += cloud_.constProps(typeIds_[iD]).mass()*(nParcels_[iD][cell]/rhoNMean_[cell]);
                        Cp += (5.0 + cloud_.constProps(typeIds_[iD]).rotationalDoF())*(nParcels_[iD][cell]/rhoNMean_[cell]);
                        Cv += (3.0 + cloud_.constProps(typeIds_[iD]).rotationalDoF())*(nParcels_[iD][cell]/rhoNMean_[cell]);
                    }
                }

                Cv_p = Cv/physicoChemical::NA.value();

                scalar gasConstant = 0.0;
                scalar gamma = 0.0;
                scalar speedOfSound = 0.0;

                if (molecularMass > VSMALL)
                {
                    gasConstant = physicoChemical::k.value()/molecularMass;
                }

                if (rhoNMean_[cell] > VSMALL)
                {
                    gamma = Cp/Cv; // gamma = cP/cV
                }

                if (rhoNMean_[cell] > VSMALL && translationalT_[cell] > VSMALL)
                {
                    speedOfSound = sqrt(gamma*gasConstant*translationalT_[cell]);
                    Ma_[cell] = mag(UMean_[cell])/speedOfSound;
                }
                else
                {
                    Ma_[cell] = 0.0;
                }

                if (measureMeanFreePath_)
                {
                    forAll(mfp_, iD)
                    {
                        label qspec = 0;

                        for (qspec=0; qspec<typeIds_.size(); ++qspec)
                        {
                            scalar dPQ = 0.5*(cloud_.constProps(typeIds_[iD]).d() + cloud_.constProps(typeIds_[qspec]).d());

                            scalar omegaPQ = 0.5*(cloud_.constProps(typeIds_[iD]).omega() + cloud_.constProps(typeIds_[qspec]).omega());

                            scalar massRatio = cloud_.constProps(typeIds_[iD]).mass()/cloud_.constProps(typeIds_[qspec]).mass();

                            if (nParcels_[qspec][cell] > VSMALL && translationalT_[cell] > VSMALL)
                            {
                                scalar nDensQ = (nParcelsXnParticle_[qspec][cell])/(mesh_.cellVolumes()[cell]*timeAvCounter_);

                                scalar reducedMass = 
                                    cloud_.constProps(typeIds_[iD]).mass()*cloud_.constProps(typeIds_[qspec]).mass()
                                    /(cloud_.constProps(typeIds_[iD]).mass() + cloud_.constProps(typeIds_[qspec]).mass());

                                //Bird, eq (4.76)
                                mfp_[iD][cell] += pi*dPQ*dPQ*nDensQ*pow(Tref_/translationalT_[cell],omegaPQ - 0.5)*sqrt(1.0 + massRatio); 

                                // //Bird, eq (4.74)
                                mcr_[iD][cell] +=
                                    (2.0*sqrt(pi)*dPQ*dPQ*nDensQ*pow(translationalT_[cell]/Tref_,1.0 - omegaPQ)*sqrt(2.0*physicoChemical::k.value()*Tref_/reducedMass)); 
                            }
                        }

                        if (mfp_[iD][cell] > VSMALL)
                        {
                            mfp_[iD][cell] = 1.0/mfp_[iD][cell];
                        }
                    }

                    MFP_[cell] = 0.0;
                    dxMFP_[cell] = 0.0;
                    MCR_[cell] = 0.0;
                    MCT_[cell] = 0.0;
                    dtMCT_[cell] = 0.0;

                    forAll(mfp_, iD)
                    {
                        if (rhoN_[cell] > VSMALL)
                        {
                            scalar nDensP = nParcelsXnParticle_[iD][cell]/(mesh_.cellVolumes()[cell]*timeAvCounter_);

                            //Bird, eq (4.77)
                            MFP_[cell] += mfp_[iD][cell]*nDensP/rhoN_[cell];

                            //Bird, eq (1.38)
                            MCR_[cell] += mcr_[iD][cell]*nDensP/rhoN_[cell];
                        }
                    }

                    if (MFP_[cell] < VSMALL)
                    {
                        MFP_[cell] = GREAT;
                    }

                    if (MCR_[cell] > VSMALL)
                    {
                        MCT_[cell] = 1.0/MCR_[cell];
                        dtMCT_[cell] = deltaT/MCT_[cell];
                    }
                    else
                    {
                        MCT_[cell] = GREAT;
                        dtMCT_[cell] = GREAT;
                    }

                    forAll(mfp_, iD)
                    {
                        mfp_[iD][cell] = 0.0;
                        mcr_[iD][cell] = 0.0;
                    }

                    if (MFP_[cell] > VSMALL)
                    {
 
                        scalar largestCellDimension = 0.0;

                        point minPoint = vector(GREAT, GREAT, GREAT);
                        point maxPoint = vector(-GREAT, -GREAT, -GREAT);
                        const List<label>& cellNodes = mesh_.cellPoints()[cell];

                        forAll(cellNodes, node) 
                        {
                            const point& cellPoint = mesh_.points()[cellNodes[node]];
                            minPoint.x() = min(minPoint.x(),cellPoint.x());
                            minPoint.y() = min(minPoint.y(),cellPoint.y());
                            minPoint.z() = min(minPoint.z(),cellPoint.z());
                            maxPoint.x() = max(maxPoint.x(),cellPoint.x());
                            maxPoint.y() = max(maxPoint.y(),cellPoint.y());
                            maxPoint.z() = max(maxPoint.z(),cellPoint.z());                
                        }

                        const boolVector& solutionDimensions = cloud_.solutionDimensions();
                        forAll(solutionDimensions, dim)
                        {
                            scalar cellDimension = (maxPoint[dim]-minPoint[dim])/cloud_.subCellLevels()[cell][dim];
                            if (solutionDimensions[dim] && largestCellDimension < cellDimension)
                            {
                                largestCellDimension = cellDimension;
                            }
                        }

                        dxMFP_[cell] = largestCellDimension/MFP_[cell];

                    }
                    else
                    {
                        dxMFP_[cell] = GREAT;
                    }
                }

                if (measureErrors_)
                {
                    if (uniGasRhoNMean_[cell] > VSMALL &&
                        Ma_[cell] > VSMALL && gamma > VSMALL
                        && Cv_p > VSMALL)
                    {
                        densityError_[cell] = 1.0/sqrt(uniGasRhoNMean_[cell]*nAvTimeSteps_);
                        velocityError_[cell] = (1.0/sqrt(uniGasRhoNMean_[cell]*nAvTimeSteps_))*(1.0/(Ma_[cell]*sqrt(gamma)));
                        temperatureError_[cell] = (1.0/sqrt(uniGasRhoNMean_[cell]*nAvTimeSteps_))*sqrt(physicoChemical::k.value()/Cv_p);
                        pressureError_[cell] = sqrt(gamma)/sqrt(uniGasRhoNMean_[cell]*nAvTimeSteps_);
                    }
                }
            }

            List<scalarField> vibTBF(mesh_.boundaryMesh().size());

            // computing boundary measurements
            forAll(boundaryCells_, j)
            {

                const polyPatch& patch = mesh_.boundaryMesh()[j];

                const labelList& bCs = boundaryCells_[j];
                
                vibTBF[j].setSize(patch.size(), 0.0);

                if (isA<wallPolyPatch>(patch))
                {
                    forAll(rhoN_.boundaryFieldRef()[j], k)
                    {
                        if (rhoNBF_[j][k] > VSMALL)
                        {

                            const vector fC = mesh_.faceCentres()[mesh_.boundaryMesh()[j].start()+k];
                            const scalar CWF = cloud_.cellWF(bCs[k]);
                            const scalar RWF = cloud_.axiRWF(fC);
                            const scalar nParticle = cloud_.nParticle()*CWF*RWF;

                            rhoN_.boundaryFieldRef()[j][k] = rhoNBF_[j][k]*nParticle/timeAvCounter_;
                            rhoM_.boundaryFieldRef()[j][k] = rhoMBF_[j][k]*nParticle/timeAvCounter_;
                            UMean_.boundaryFieldRef()[j][k] = momentumBF_[j][k]*nParticle/(rhoM_.boundaryFieldRef()[j][k]*timeAvCounter_);

                            scalar rhoMMean = rhoMBF_[j][k]*nParticle/timeAvCounter_;
                            scalar linearKEMean = linearKEBF_[j][k]*nParticle/timeAvCounter_;
                            scalar rhoNMean = rhoNBF_[j][k]*nParticle/timeAvCounter_;
                            translationalT_.boundaryFieldRef()[j][k] =
                                2.0/(3.0*physicoChemical::k.value()*rhoNMean)
                                *(linearKEMean - 0.5*rhoMMean*(UMean_.boundaryFieldRef()[j][k] & UMean_.boundaryFieldRef()[j][k]));

                        }
                        else
                        {
                            rhoN_.boundaryFieldRef()[j][k] = 0.0;
                            rhoM_.boundaryFieldRef()[j][k] = 0.0;
                            translationalT_.boundaryFieldRef()[j][k] = 0.0;
                            UMean_.boundaryFieldRef()[j][k] = vector::zero;
                        }

                        if (rotationalDofBF_[j][k] > VSMALL)
                        {
                            rotationalT_.boundaryFieldRef()[j][k] =
                                (2.0/physicoChemical::k.value())
                               *(rotationalEBF_[j][k]/rotationalDofBF_[j][k]);
                        }
                        else
                        {
                            rotationalT_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        // electronic temperature
                        scalar totalEDof = 0.0;
                        scalar elecT = 0.0;

                        electronicT_.boundaryFieldRef()[j][k] = elecT;

                        scalar nRotDof = 0.0;

                        if (rhoNBF_[j][k] > VSMALL)
                        {
                            nRotDof = rotationalDofBF_[j][k]/rhoNBF_[j][k];
                        }

                        overallT_.boundaryFieldRef()[j][k] =
                            (
                                (3.0*translationalT_.boundaryFieldRef()[j][k])
                              + (nRotDof*rotationalT_.boundaryFieldRef()[j][k])
                              + (
                                    totalvDofBF_[j][k]
                                   *vibrationalT_.boundaryFieldRef()[j][k]
                                )
                              + (totalEDof*elecT)
                            )
                           /(3.0 + nRotDof + totalvDofBF_[j][k] + totalEDof);

                        totalvDofBF_[j][k] = 0.0;

                        q_.boundaryFieldRef()[j][k] = qBF_[j][k]/timeAvCounter_;
                        fD_.boundaryFieldRef()[j][k] = fDBF_[j][k]/timeAvCounter_;

                    }

                    p_.boundaryFieldRef()[j] = fD_.boundaryFieldRef()[j] & n_[j];

                    tau_.boundaryFieldRef()[j] =
                        sqrt
                        (
                            sqr(fD_.boundaryFieldRef()[j] & t1_[j])
                          + sqr(fD_.boundaryFieldRef()[j] & t2_[j])
                        );
                }
            }

            forAll(boundaryCells_, j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];

                const labelList& bCs = boundaryCells_[j];

                forAll(bCs, k)
                {
                    if
                    (
                        isA<polyPatch>(patch)
                     && !isA<emptyPolyPatch>(patch)
                     && !isA<cyclicPolyPatch>(patch)
                    )
                    {

                        uniGasRhoNMean_.boundaryFieldRef()[j][k] = uniGasRhoNMean_[bCs[k]];
                        rhoN_.boundaryFieldRef()[j][k] = rhoN_[bCs[k]];
                        rhoM_.boundaryFieldRef()[j][k] = rhoM_[bCs[k]];
                        Ma_.boundaryFieldRef()[j][k] = Ma_[bCs[k]];

                        if (measureMeanFreePath_)
                        {
                            MFP_.boundaryFieldRef()[j][k] = MFP_[bCs[k]];
                            dxMFP_.boundaryFieldRef()[j][k] = dxMFP_[bCs[k]];
                            MCR_.boundaryFieldRef()[j][k] = MCR_[bCs[k]];
                            MCT_.boundaryFieldRef()[j][k] = MCT_[bCs[k]];
                            dtMCT_.boundaryFieldRef()[j][k] = dtMCT_[bCs[k]];
                        }

                        if (!isA<wallPolyPatch>(patch))
                        {
                            translationalT_.boundaryFieldRef()[j][k] = translationalT_[bCs[k]];
                            rotationalT_.boundaryFieldRef()[j][k] = rotationalT_[bCs[k]];
                            vibrationalT_.boundaryFieldRef()[j][k] = vibrationalT_[bCs[k]];
                            overallT_.boundaryFieldRef()[j][k] = overallT_[bCs[k]];
                            p_.boundaryFieldRef()[j][k] = p_[bCs[k]];
                            UMean_.boundaryFieldRef()[j][k] = UMean_[bCs[k]];
                        }
                    }
                }
            }

            if (measureMeanFreePath_)
            {
                MFP_.write();
                dxMFP_.write();
                MCR_.write();
                MCT_.write();
                dtMCT_.write();
            }


            if (measureErrors_)
            {
                densityError_.write();
                velocityError_.write();
                temperatureError_.write();
                pressureError_.write();
            }

            uniGasRhoNMean_.write();
            rhoN_.write();
            rhoM_.write();
            p_.write();
            translationalT_.write();
            rotationalT_.write();
            vibrationalT_.write();
            electronicT_.write();
            overallT_.write();
            q_.write();
            tau_.write();
            Ma_.write();
            UMean_.write();
            fD_.write();

        }

        // Reset
        if (resetFieldsAtOutput_ && mesh_.time().value() < resetFieldsAtOutputUntilTime_+0.5*deltaT)
        {

            nAvTimeSteps_ = 0;

            timeAvCounter_ = 0.0;

            forAll(rhoNMean_, cell)
            {
                rhoNMean_[cell] = 0.0;
                rhoMMean_[cell] = 0.0;
                linearKEMean_[cell] = 0.0;
                momentumMean_[cell] = vector::zero;
                rotationalEMean_[cell] = 0.0;
                rotationalDofMean_[cell] = 0.0;
                rhoNMeanInt_[cell] = 0.0;
                molsElec_[cell] = 0.0,
                rhoNMeanXnParticle_[cell] = 0.0;
                rhoMMeanXnParticle_[cell] = 0.0;
                momentumMeanXnParticle_[cell] = vector::zero;
                linearKEMeanXnParticle_[cell] = 0.0;
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

            // reset boundary information
            forAll(rhoNBF_, j)
            {
                rhoNBF_[j] = 0.0;
                rhoMBF_[j] = 0.0;
                linearKEBF_[j] = 0.0;
                speciesRhoNIntBF_[j] = 0.0;
                speciesRhoNElecBF_[j] = 0.0;
                rotationalEBF_[j] = 0.0;
                rotationalDofBF_[j] = 0.0;
                totalvDofBF_[j] = 0.0;
                qBF_[j] = 0.0;
                fDBF_[j] = vector::zero;
                momentumBF_[j] = vector::zero;
            }

            forAll(speciesRhoNBF_, i)
            {
                forAll(speciesRhoNBF_[i], j)
                {
                    speciesRhoNBF_[i][j] = 0.0;
                    vibrationalEBF_[i][j] = 0.0;
                    electronicEBF_[i][j] = 0.0;
                    mccSpeciesBF_[i][j] = 0.0;
                }
            }
        }

    }

    if (averagingAcrossManyRuns_)
    {
        writeOut();
    }

}


void Foam::uniGasVolFields::writeField()
{}


void Foam::uniGasVolFields::updateProperties(const dictionary& dict)
{
    
    // The main properties should be updated first
    uniGasField::updateProperties(dict);
    
    propsDict_ = dict.subDict(typeName + "Properties");

    propsDict_.readIfPresent("Tref", Tref_);

    propsDict_.readIfPresent("measureErrors", measureErrors_);

    propsDict_.readIfPresent("densityOnly", densityOnly_);

    propsDict_.readIfPresent("measureMeanFreePath", measureMeanFreePath_);

    propsDict_.readIfPresent("averagingAcrossManyRuns", averagingAcrossManyRuns_);

}
// ************************************************************************* //