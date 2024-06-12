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

\*---------------------------------------------------------------------------*/

#include "uniGasGeneralBoundary.H"
#include "fvMesh.H"
#include "graph.H"
#include "mathematicalConstants.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uniGasGeneralBoundary, 0);
defineRunTimeSelectionTable(uniGasGeneralBoundary, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasGeneralBoundary::uniGasGeneralBoundary
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
:
    uniGasBoundaryBase(mesh, cloud, dict, "general"),
    faces_(),
    patchSurfaceArea_(0.0),
    cells_(),
    accumulatedParcelsToInsert_()
{
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];

    // Initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    // Loop through all faces and set the boundary cells
    // - no conflict with parallelisation because the faces are unique

    for (label i = 0; i < patch.size(); ++i)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }

    if (Pstream::parRun())
    {
        reduce(patchSurfaceArea_, sumOp<scalar>());
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::uniGasGeneralBoundary> Foam::uniGasGeneralBoundary::New
(
    const polyMesh& mesh,
    uniGasCloud& cloud,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("boundaryModel"));

    Info<< "Selecting uniGasGeneralBoundaryModel " << modelType << nl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "boundaryModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<uniGasGeneralBoundary>(cstrIter()(mesh, cloud, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniGasGeneralBoundary::computeParcelsToInsert
(
    const scalarList& numDen,
    const scalar& transT,
    const vector& velocity
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    scalar sqrtPi = sqrt(pi);

    // compute parcels to insert
    forAll(accumulatedParcelsToInsert_, i)
    {
        const label typeId = typeIds_[i];
        const scalar& mass = cloud_.constProps(typeId).mass();

        forAll(accumulatedParcelsToInsert_[i], f)
        {
            const label faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT,
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity & -sF/fA)/mostProbableSpeed;

            // From Bird eqn 4.22
            scalar CWF = cloud_.cellWF(cells_[f]);
            scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);

            accumulatedParcelsToInsert_[i][f] +=
            (
                fA*numDen[i]*deltaT*mostProbableSpeed
               *(
                    exp(-sqr(sCosTheta))
                  + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                )
            )
            /(2.0*sqrtPi*cloud_.nParticle()*CWF*RWF);

        }
    }
}

void Foam::uniGasGeneralBoundary::computeParcelsToInsert
(
    const scalarList& numDen,
    const scalar& transT,
    const vector& velocity,
    const vector& heatFlux,
    const tensor& stress
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    scalar sqrtPi = sqrt(pi);

    scalar pressure=0.0;
    forAll(accumulatedParcelsToInsert_, i)
    {
        pressure+=numDen[i];
    }    
    pressure*=physicoChemical::k.value()*transT;

    // compute parcels to insert
    forAll(accumulatedParcelsToInsert_, i)
    {
        const label typeId = typeIds_[i];
        const scalar& mass = cloud_.constProps(typeId).mass();

        forAll(accumulatedParcelsToInsert_[i], f)
        {
            const label faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT,
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity & -sF/fA)/mostProbableSpeed;

            //  compute heatflux and stress normal to the boundary face
            scalar heatFluxNormal = ( heatFlux & -sF/fA );

            scalar stressNormal = ((stress & -sF/fA) & -sF/fA);

            // From Bird eqn 4.22
            scalar CWF = cloud_.cellWF(cells_[f]);
            scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);

            accumulatedParcelsToInsert_[i][f] +=
            (
                fA*numDen[i]*deltaT*mostProbableSpeed
               *(
                    exp(-sqr(sCosTheta))
                  *(1.0-0.5*stressNormal/pressure-0.4*heatFluxNormal*sCosTheta/pressure/mostProbableSpeed)
                  + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                )
            )
            /(2.0*sqrtPi*cloud_.nParticle()*CWF*RWF);
        }
    }
}

void Foam::uniGasGeneralBoundary::computeParcelsToInsert
(
    const List<scalarField>& numDen,
    const scalarField& transT,
    const vectorField& velocity
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    const scalar sqrtPi = sqrt(pi);

    // Compute parcels to insert
    forAll(accumulatedParcelsToInsert_, i)
    {
        const label typeId = typeIds_[i];
        const scalar& mass = cloud_.constProps(typeId).mass();

        forAll(accumulatedParcelsToInsert_[i], f)
        {
            const label faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT[f],
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity[f] & -sF/fA)
                /mostProbableSpeed;

            // From Bird eqn 4.22
            scalar CWF = cloud_.cellWF(cells_[f]);
            scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);

            accumulatedParcelsToInsert_[i][f] +=
            (
                fA*numDen[i][f]*deltaT*mostProbableSpeed
               *(
                    exp(-sqr(sCosTheta))
                  + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                )
            )
            /(2.0*sqrtPi*cloud_.nParticle()*CWF*RWF);

        }
    }
}

void Foam::uniGasGeneralBoundary::computeParcelsToInsert
(
    const List<scalarField>& numDen,
    const scalarField& transT,
    const vectorField& velocity,
    const vectorField& heatFlux,
    const tensorField& stress
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    const scalar sqrtPi = sqrt(pi);

    // Compute parcels to insert
    forAll(accumulatedParcelsToInsert_, i)
    {
        const label typeId = typeIds_[i];
        const scalar& mass = cloud_.constProps(typeId).mass();

        forAll(accumulatedParcelsToInsert_[i], f)
        {
            const label faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            scalar pressure=0.0;
            forAll(accumulatedParcelsToInsert_, i)
            {
                pressure+=numDen[i][f];
            }    
            pressure*=physicoChemical::k.value()*transT[f];

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT[f],
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity[f] & -sF/fA)
                /mostProbableSpeed;

            //  compute heatflux and stress normal to the boundary face
            scalar heatFluxNormal = ( heatFlux[f] & -sF/fA );

            scalar stressNormal = ((stress[f] & -sF/fA) & -sF/fA);

            // From Bird eqn 4.22
            scalar CWF = cloud_.cellWF(cells_[f]);
            scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);

            accumulatedParcelsToInsert_[i][f] +=
            (
                fA*numDen[i][f]*deltaT*mostProbableSpeed
               *(
                    exp(-sqr(sCosTheta))
                  *(1.0-0.5*stressNormal/pressure-0.4*heatFluxNormal*sCosTheta/pressure/mostProbableSpeed)
                  + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                )
            )
            /(2.0*sqrtPi*cloud_.nParticle()*CWF*RWF);
        }
    }
}

void Foam::uniGasGeneralBoundary::computeParcelsToInsert
(
    const scalar& numDen,
    const scalarField& molFractions,
    const scalar& transT,
    const vectorField& velocity
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    scalar sqrtPi = sqrt(pi);
    
    forAll(accumulatedParcelsToInsert_, iD)
    {
        const label typeId = typeIds_[iD];

        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const label faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar& fA = mag(sF);

            const scalar& pMass = cloud_.constProps(typeId).mass();

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT,
                    pMass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity[f] & -sF/fA)
                                    /mostProbableSpeed;

            // From Bird eqn 4.22
            scalar CWF = cloud_.cellWF(cells_[f]);
            scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faceI]);

            accumulatedParcelsToInsert_[iD][f] +=
                molFractions[iD]
               *(
                    fA*numDen*deltaT*mostProbableSpeed
                   *(
                        exp(-sqr(sCosTheta))
                      + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi*cloud_.nParticle()*CWF*RWF);
        }
    }
}

void Foam::uniGasGeneralBoundary::computeParcelsToInsert
(
    const scalarField& numDen,
    const scalarField& molFractions,
    const scalarField& transT,
    const vectorField& velocity
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    scalar sqrtPi = sqrt(pi);
    
    forAll(accumulatedParcelsToInsert_, iD)
    {
        const label typeId = typeIds_[iD];

        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const vector& sF = mesh_.faceAreas()[faces_[f]];
            const scalar fA = mag(sF);

            scalar mass = cloud_.constProps(typeId).mass();

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT[f],
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity[f] & -sF/fA)/mostProbableSpeed;

            // From Bird eqn 4.22
            scalar CWF = cloud_.cellWF(cells_[f]);
            scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faces_[f]]);

            accumulatedParcelsToInsert_[iD][f] +=
                molFractions[iD]
               *(
                    fA*numDen[f]*deltaT*mostProbableSpeed
                   *(
                        exp(-sqr(sCosTheta))
                      + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
               /(2.0*sqrtPi*cloud_.nParticle()*CWF*RWF);
        }
    }
}

void Foam::uniGasGeneralBoundary::computeParcelsToInsert
(
    const List<scalarField>& numDen,
    const scalar& transT,
    const vectorField& velocity
)
{
    const scalar& deltaT = mesh_.time().deltaTValue();
    scalar sqrtPi = sqrt(pi);
    
    forAll(accumulatedParcelsToInsert_, iD)
    {
        const label typeId = typeIds_[iD];

        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const vector& sF = mesh_.faceAreas()[faces_[f]];
            const scalar fA = mag(sF);

            scalar mass = cloud_.constProps(typeId).mass();

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    transT,
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            scalar sCosTheta = (velocity[f] & -sF/fA )/mostProbableSpeed;

            scalar CWF = cloud_.cellWF(cells_[f]);
            scalar RWF = cloud_.axiRWF(cloud_.mesh().faceCentres()[faces_[f]]);

            accumulatedParcelsToInsert_[iD][f] +=
            (
                fA*numDen[iD][f]*deltaT*mostProbableSpeed
               *(
                    exp(-sqr(sCosTheta))
                  + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                )
            )
           /(2.0*sqrtPi*cloud_.nParticle()*CWF*RWF);
        }
    } 
}


void Foam::uniGasGeneralBoundary::insertParcels
(
    const scalar& transT,
    const scalar& rotT,
    const scalar& vibT,
    const scalar& elecT,
    const vector& velocity
)
{

    Random& rndGen = cloud_.rndGen();

    // insert pacels
    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const label cellI = cells_[f];
        const vector& fC = mesh_.faceCentres()[faceI];
        const vector& sF = mesh_.faceAreas()[faces_[f]];
        scalar fA = mag(sF);

        List<tetIndices> faceTets =
            polyMeshTetDecomposition::faceTetIndices
            (
                mesh_,
                faceI,
                cellI
            );

        // Cumulative triangle area fractions
        List<scalar> cTriAFracs(faceTets.size(), Zero);

        cTriAFracs[0] = faceTets[0].faceTri(mesh_).mag()/fA;
        for (label trii = 1; trii < cTriAFracs.size(); ++trii)
        {
            cTriAFracs[trii] =
                cTriAFracs[trii-1] + faceTets[trii].faceTri(mesh_).mag()/fA;
        }

        // Force the last area fraction value to 1.0 to avoid any
        // rounding/non - flat face errors giving a value < 1.0
        cTriAFracs.last() = 1.0;

        // Normal unit vector *negative* so normal is pointing into the
        // domain
        vector n = sF;
        n /= -mag(n);

        // Wall tangential unit vector. Use the direction between the
        // face centre and the first vertex in the list
        vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
        t1 /= mag(t1);

        // Other tangential unit vector.  Rescaling in case face is not
        // flat and n and t1 aren't perfectly orthogonal
        vector t2 = n^t1;
        t2 /= mag(t2);

        forAll(typeIds_, iD)
        {
            const label typeId = typeIds_[iD];

            scalar& faceAccumulator = accumulatedParcelsToInsert_[iD][f];

            // Number of whole particles to insert
            label nParcelsToInsert = max(label(faceAccumulator), 0);

            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of faceAccumulator
            if ((faceAccumulator - nParcelsToInsert) > rndGen.sample01<scalar>())
            {
                ++nParcelsToInsert;
            }

            faceAccumulator -= nParcelsToInsert;

            const scalar mass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nParcelsToInsert; ++i)
            {
                // Choose a triangle to insert on, based on their relative
                // area

                scalar triSelection = rndGen.sample01<scalar>();

                // Selected triangle
                label selectedTriI = -1;

                forAll(cTriAFracs, triI)
                {
                    selectedTriI = triI;

                    if (cTriAFracs[triI] >= triSelection)
                    {
                        break;
                    }
                }

                // Randomly distribute the points on the triangle.

                const tetIndices& faceTetIs = faceTets[selectedTriI];

                point p = faceTetIs.faceTri(mesh_).randomPoint(rndGen);

                // Velocity generation
                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        transT,
                        mass
                    )
                );

                scalar sCosTheta = (velocity & n)/mostProbableSpeed;

                // Coefficients required for Bird eqn 12.5
                scalar uNormProbCoeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                scalar uNormProbCoeffB =
                    0.5
                   *(
                        1.0
                      + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                    );

                // Equivalent to the QA value in Bird's DSMC3.FOR
                scalar randomScaling = 3.0;

                if (sCosTheta < -3)
                {
                    randomScaling = mag(sCosTheta) + 1;
                }

                scalar P = -1;

                // Normalised candidates for the normal direction velocity
                // component
                scalar uNormal;
                scalar uNormalThermal;

                if (mag(velocity & n) > VSMALL)
                {
                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P =
                                2.0*uNormal/uNormProbCoeffA
                               *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.sample01<scalar>());
                }
                else
                {
                    uNormal = sqrt(-log(rndGen.sample01<scalar>()));
                }

                vector U =
                    sqrt(physicoChemical::k.value()
                   *transT/mass)
                   *(
                        rndGen.GaussNormal<scalar>()*t1
                        + rndGen.GaussNormal<scalar>()*t2
                    )
                  + (t1 & velocity)*t1
                  + (t2 & velocity)*t2
                  + mostProbableSpeed*uNormal*n;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    rotT,
                    cloud_.constProps(typeId).rotationalDoF()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    vibT,
                    cloud_.constProps(typeId).vibrationalDoF(),
                    typeId
                );

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    elecT,
                    cloud_.constProps(typeId).degeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList(),
                    typeId
                );

                label newParcel = 1;

                scalar CWF = cloud_.cellWF(cellI);
                scalar RWF = cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);
                
                // Apply tracking correction towards cell centre
                p += VSMALL*(mesh_.cellCentres()[cellI] - p);

                cloud_.addNewParcel
                (
                    p,
                    U,
                    CWF,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    typeId,
                    newParcel,
                    vibLevel
                );

            }
        }
    }
}

void Foam::uniGasGeneralBoundary::insertParcels
(
    const scalarList& numDen,
    const scalar& transT,
    const scalar& rotT,
    const scalar& vibT,
    const scalar& elecT,
    const vector& velocity,
    const vector& heatFlux,
    const tensor& stress
)
{
    Random& rndGen = cloud_.rndGen();

    // insert pacels
    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const label cellI = cells_[f];
        const vector& fC = mesh_.faceCentres()[faceI];
        const vector& sF = mesh_.faceAreas()[faces_[f]];
        scalar fA = mag(sF);

        List<tetIndices> faceTets =
            polyMeshTetDecomposition::faceTetIndices
            (
                mesh_,
                faceI,
                cellI
            );

        // Cumulative triangle area fractions
        List<scalar> cTriAFracs(faceTets.size(), Zero);

        cTriAFracs[0] = faceTets[0].faceTri(mesh_).mag()/fA;
        for (label trii = 1; trii < cTriAFracs.size(); ++trii)
        {
            cTriAFracs[trii] =
                cTriAFracs[trii-1] + faceTets[trii].faceTri(mesh_).mag()/fA;
        }

        // Force the last area fraction value to 1.0 to avoid any
        // rounding/non - flat face errors giving a value < 1.0
        cTriAFracs.last() = 1.0;

        // Normal unit vector *negative* so normal is pointing into the
        // domain
        vector n = sF;
        n /= -mag(n);

        // Wall tangential unit vector. Use the direction between the
        // face centre and the first vertex in the list
        vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
        t1 /= mag(t1);

        // Other tangential unit vector.  Rescaling in case face is not
        // flat and n and t1 aren't perfectly orthogonal
        vector t2 = n^t1;
        t2 /= mag(t2);

        //Compute pressure of mixture based on number densities of all species
        scalar pressure=0.0;
        forAll(typeIds_, iD)
        {
            pressure+=numDen[iD];
        }
        pressure*=physicoChemical::k.value()*transT;

        forAll(typeIds_, iD)
        {
            const label typeId = typeIds_[iD];

            scalar& faceAccumulator = accumulatedParcelsToInsert_[iD][f];

            // Number of whole particles to insert
            label nParcelsToInsert = max(label(faceAccumulator), 0);

            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of faceAccumulator
            if ((faceAccumulator - nParcelsToInsert) > rndGen.sample01<scalar>())
            {
                ++nParcelsToInsert;
            }

            faceAccumulator -= nParcelsToInsert;

            const scalar mass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nParcelsToInsert; ++i)
            {
                // Choose a triangle to insert on, based on their relative
                // area

                scalar triSelection = rndGen.sample01<scalar>();

                // Selected triangle
                label selectedTriI = -1;

                forAll(cTriAFracs, triI)
                {
                    selectedTriI = triI;

                    if (cTriAFracs[triI] >= triSelection)
                    {
                        break;
                    }
                }

                // Randomly distribute the points on the triangle.

                const tetIndices& faceTetIs = faceTets[selectedTriI];

                point p = faceTetIs.faceTri(mesh_).randomPoint(rndGen);

                // Velocity generation
               scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        transT,
                        mass
                    )
                );

                // Compute breakdown parameter
                scalar maxHeatFlux=-1.0;
                forAll(heatFlux,i) 
                {
                    if (maxHeatFlux < fabs(heatFlux[i])) 
                    {
                        maxHeatFlux = fabs(heatFlux[i]);
                    }
                }
                maxHeatFlux = 2.0*maxHeatFlux/(pressure*mostProbableSpeed);

                scalar maxStress=-1.0;
                forAll(stress,i) 
                {
                    if (maxStress < fabs(stress[i])) 
                    {
                        maxStress = fabs(stress[i]);
                    }
                }
                maxStress = maxStress/pressure;

                scalar breakdownParameter = max(maxHeatFlux,maxStress);

                // Compute amplitude parameter
                scalar amplitudeParameter = 1.0 + 60.0*breakdownParameter;
                scalar gamma=0.0;

                scalar sCosTheta = (velocity & n)/mostProbableSpeed;

                scalar uNormalCoeffA=4.0;
                scalar uNormalCoeffB=5.0;
                scalar lowerBound=min(sCosTheta-uNormalCoeffA,-uNormalCoeffB);
                scalar upperBound=min(sCosTheta,uNormalCoeffB);
                scalar uNormalMax=0.5*(sCosTheta-sqrt(sqr(sCosTheta)+2.0));

                // Select a velocity using rejection algorithm in [A. L. Garcia, 2006]
                vector U=vector::zero;
                do 
                {

                    scalar uNormal;
                    
                    if(abs(velocity & n) > VSMALL)
                    {
                        // Select a velocity using box envelope method described in [A. L. Garcia, 1998]
                        do
                        {
                            
                            uNormal=lowerBound+rndGen.sample01<scalar>()*(upperBound-lowerBound);

                        } while ( (sCosTheta-uNormal)/(sCosTheta-uNormalMax)*exp(sqr(uNormalMax)-sqr(uNormal)) < rndGen.sample01<scalar>());
                    }
                    else
                    {
                        uNormal = -sqrt(-log(rndGen.sample01<scalar>()));
                    }

                    U = rndGen.GaussNormal<scalar>()/sqrt(2.0)*t1
                        + rndGen.GaussNormal<scalar>()/sqrt(2.0)*t2
                        - uNormal*n;

                    gamma=1.0+(2.0/mostProbableSpeed*(heatFlux.x()*U.x()+heatFlux.y()*U.y()+heatFlux.z()*U.z())*(2.0/5.0*(sqr(U.x())+sqr(U.y())+sqr(U.z()))-1.0)
                        -2.0*(stress.xy()*U.x()*U.y()+stress.xz()*U.x()*U.z()+stress.yz()*U.y()*U.z())
                        -stress.xx()*(sqr(U.x())-sqr(U.z()))-stress.yy()*(sqr(U.y())-sqr(U.z())))/pressure;

                } while(amplitudeParameter*rndGen.sample01<scalar>()>gamma);

                U=mostProbableSpeed*U+velocity;
                /////////////////////////////////////////////////////////////////////////

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    rotT,
                    cloud_.constProps(typeId).rotationalDoF()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    vibT,
                    cloud_.constProps(typeId).vibrationalDoF(),
                    typeId
                );

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    elecT,
                    cloud_.constProps(typeId).degeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList(),
                    typeId
                );

                label newParcel = 1;

                scalar CWF = cloud_.cellWF(cellI);
                scalar RWF = cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);
                
                // Apply tracking correction towards cell centre
                p += VSMALL*(mesh_.cellCentres()[cellI] - p);

                cloud_.addNewParcel
                (
                    p,
                    U,
                    CWF,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    typeId,
                    newParcel,
                    vibLevel
                );

            }
        }
    }
}

void Foam::uniGasGeneralBoundary::insertParcels
(
    const scalarField& transT,
    const scalarField& rotT,
    const scalarField& vibT,
    const scalarField& elecT,
    const vectorField& velocity
)
{
    Random& rndGen = cloud_.rndGen();

    // Insert parcels
    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const label cellI = cells_[f];
        const vector& fC = mesh_.faceCentres()[faceI];
        const vector& sF = mesh_.faceAreas()[faces_[f]];
        scalar fA = mag(sF);
        const vector& faceVelocity = velocity[f];
        const scalar faceTranslationalTemperature = transT[f];
        const scalar faceRotationalTemperature = rotT[f];
        const scalar faceVibrationalTemperature = vibT[f];
        const scalar faceElectronicTemperature = elecT[f];

        List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
        (
            mesh_,
            faceI,
            cellI
        );

        // Cumulative triangle area fractions
        List<scalar> cTriAFracs(faceTets.size(), Zero);

        cTriAFracs[0] = faceTets[0].faceTri(mesh_).mag()/fA;
        for (label trii = 1; trii < cTriAFracs.size(); ++trii)
        {
            cTriAFracs[trii] =
                cTriAFracs[trii-1] + faceTets[trii].faceTri(mesh_).mag()/fA;
        }

        // Force the last area fraction value to 1.0 to avoid any
        // rounding/non - flat face errors giving a value < 1.0
        cTriAFracs.last() = 1.0;

        // Normal unit vector *negative* so normal is pointing into the
        // domain
        vector n = sF;
        n /= -mag(n);

        // Wall tangential unit vector. Use the direction between the
        // face centre and the first vertex in the list
        vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
        t1 /= mag(t1);

        // Other tangential unit vector.  Rescaling in case face is not
        // flat and n and t1 aren't perfectly orthogonal
        vector t2 = n^t1;
        t2 /= mag(t2);

        forAll(typeIds_, iD)
        {
            const label typeId = typeIds_[iD];

            scalar& faceAccumulator = accumulatedParcelsToInsert_[iD][f];

            // Number of whole particles to insert
            label nParcelsToInsert = max(label(faceAccumulator), 0);

            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of faceAccumulator
            if ((faceAccumulator - nParcelsToInsert) > rndGen.sample01<scalar>())
            {
                ++nParcelsToInsert;
            }

            faceAccumulator -= nParcelsToInsert;

            scalar mass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nParcelsToInsert; ++i)
            {
                // Choose a triangle to insert on, based on their relative
                // area

                scalar triSelection = rndGen.sample01<scalar>();

                // Selected triangle
                label selectedTriI = -1;

                forAll(cTriAFracs, triI)
                {
                    selectedTriI = triI;

                    if (cTriAFracs[triI] >= triSelection)
                    {
                        break;
                    }
                }

                // Randomly distribute the points on the triangle.

                const tetIndices& faceTetIs = faceTets[selectedTriI];

                point p = faceTetIs.faceTri(mesh_).randomPoint(rndGen);

                // Velocity generation

                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        faceTranslationalTemperature,
                        mass
                    )
                );

                scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                // Coefficients required for Bird eqn 12.5
                scalar uNormProbCoeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                scalar uNormProbCoeffB =
                    0.5
                   *(
                        1.0
                      + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                    );

                // Equivalent to the QA value in Bird's DSMC3.FOR
                scalar randomScaling = 3.0;

                if (sCosTheta < -3)
                {
                    randomScaling = mag(sCosTheta) + 1;
                }

                scalar P = -1;

                // Normalised candidates for the normal direction velocity
                // component
                scalar uNormal;
                scalar uNormalThermal;

                // Select a velocity using Bird eqn 12.5
                do
                {
                    uNormalThermal =
                        randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                    uNormal = uNormalThermal + sCosTheta;

                    if (uNormal < 0.0)
                    {
                        P = -1;
                    }
                    else
                    {
                        P =
                            2.0*uNormal/uNormProbCoeffA
                           *exp(uNormProbCoeffB - sqr(uNormalThermal));
                    }

                } while (P < rndGen.sample01<scalar>());

                vector U =
                    sqrt(physicoChemical::k.value()
                       *faceTranslationalTemperature/mass)
                       *(
                            rndGen.GaussNormal<scalar>()*t1
                          + rndGen.GaussNormal<scalar>()*t2
                        )
                      + (t1 & faceVelocity)*t1
                      + (t2 & faceVelocity)*t2
                      + mostProbableSpeed*uNormal*n;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceRotationalTemperature,
                    cloud_.constProps(typeId).rotationalDoF()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    faceVibrationalTemperature,
                    cloud_.constProps(typeId).vibrationalDoF(),
                    typeId
                );

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceElectronicTemperature,
                    cloud_.constProps(typeId).degeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList(),
                    typeId
                );


                label newParcel = 1;

                scalar CWF = cloud_.cellWF(cellI);
                scalar RWF = cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);
                
                // Apply tracking correction towards cell centre
                p += VSMALL*(mesh_.cellCentres()[cellI] - p);

                cloud_.addNewParcel
                (
                    p,
                    U,
                    CWF,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    typeId,
                    newParcel,
                    vibLevel
                );

            }
        }
    }
}

void Foam::uniGasGeneralBoundary::insertParcels
(
    const List<scalarField>& numDen,
    const scalarField& transT,
    const scalarField& rotT,
    const scalarField& vibT,
    const scalarField& elecT,
    const vectorField& velocity,
    const vectorField& heatFlux,
    const tensorField& stress
)
{
    Random& rndGen = cloud_.rndGen();

    // Insert parcels
    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const label cellI = cells_[f];
        const vector& fC = mesh_.faceCentres()[faceI];
        const vector& sF = mesh_.faceAreas()[faces_[f]];
        scalar fA = mag(sF);
        const vector& faceVelocity = velocity[f];
        const scalar faceTranslationalTemperature = transT[f];
        const scalar faceRotationalTemperature = rotT[f];
        const scalar faceVibrationalTemperature = vibT[f];
        const scalar faceElectronicTemperature = elecT[f];
        const vector faceHeatFlux = heatFlux[f];
        const tensor faceStress = stress[f];

        List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
        (
            mesh_,
            faceI,
            cellI
        );

        // Cumulative triangle area fractions
        List<scalar> cTriAFracs(faceTets.size(), Zero);

        cTriAFracs[0] = faceTets[0].faceTri(mesh_).mag()/fA;
        for (label trii = 1; trii < cTriAFracs.size(); ++trii)
        {
            cTriAFracs[trii] =
                cTriAFracs[trii-1] + faceTets[trii].faceTri(mesh_).mag()/fA;
        }

        // Force the last area fraction value to 1.0 to avoid any
        // rounding/non - flat face errors giving a value < 1.0
        cTriAFracs.last() = 1.0;

        // Normal unit vector *negative* so normal is pointing into the
        // domain
        vector n = sF;
        n /= -mag(n);

        // Wall tangential unit vector. Use the direction between the
        // face centre and the first vertex in the list
        vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
        t1 /= mag(t1);

        // Other tangential unit vector.  Rescaling in case face is not
        // flat and n and t1 aren't perfectly orthogonal
        vector t2 = n^t1;
        t2 /= mag(t2);

        //Compute pressure of mixture based on number densities of all species
        scalar pressure=0.0;
        forAll(typeIds_, iD)
        {
            pressure+=numDen[iD][f];
        }
        pressure*=physicoChemical::k.value()*faceTranslationalTemperature;

        forAll(typeIds_, iD)
        {
            const label typeId = typeIds_[iD];

            scalar& faceAccumulator = accumulatedParcelsToInsert_[iD][f];

            // Number of whole particles to insert
            label nParcelsToInsert = max(label(faceAccumulator), 0);

            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of faceAccumulator
            if ((faceAccumulator - nParcelsToInsert) > rndGen.sample01<scalar>())
            {
                ++nParcelsToInsert;
            }

            faceAccumulator -= nParcelsToInsert;

            scalar mass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nParcelsToInsert; ++i)
            {
                // Choose a triangle to insert on, based on their relative
                // area

                scalar triSelection = rndGen.sample01<scalar>();

                // Selected triangle
                label selectedTriI = -1;

                forAll(cTriAFracs, triI)
                {
                    selectedTriI = triI;

                    if (cTriAFracs[triI] >= triSelection)
                    {
                        break;
                    }
                }

                // Randomly distribute the points on the triangle.

                const tetIndices& faceTetIs = faceTets[selectedTriI];

                point p = faceTetIs.faceTri(mesh_).randomPoint(rndGen);

                // Velocity generation

                //Compute pressure of mixture based on number densities of all species

                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        faceTranslationalTemperature,
                        mass
                    )
                );

                // Compute breakdown parameter
                scalar maxHeatFlux=-1.0;
                forAll(faceHeatFlux,i) 
                {
                    if (maxHeatFlux < fabs(faceHeatFlux[i])) 
                    {
                        maxHeatFlux = fabs(faceHeatFlux[i]);
                    }
                }
                maxHeatFlux = 2.0*maxHeatFlux/(pressure*mostProbableSpeed);

                scalar maxStress=-1.0;
                forAll(faceStress,i) 
                {
                    if (maxStress < fabs(faceStress[i])) 
                    {
                        maxStress = fabs(faceStress[i]);
                    }
                }
                maxStress = maxStress/pressure;

                scalar breakdownParameter = max(maxHeatFlux,maxStress);

                // Compute amplitude parameter
                scalar amplitudeParameter = 1.0 + 60.0*breakdownParameter;
                scalar gamma=0.0;

                scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                scalar uNormalCoeffA=4.0;
                scalar uNormalCoeffB=5.0;
                scalar lowerBound=min(sCosTheta-uNormalCoeffA,-uNormalCoeffB);
                scalar upperBound=min(sCosTheta,uNormalCoeffB);
                scalar uNormalMax=0.5*(sCosTheta-sqrt(sqr(sCosTheta)+2.0));

                // Select a velocity using rejection algorithm in [A. L. Garcia, 2006]
                vector U=vector::zero;
                do 
                {

                    scalar uNormal;
                    
                    if(abs(faceVelocity & n) > VSMALL)
                    {
                        // Select a velocity using box envelope method described in [A. L. Garcia, 1998]
                        do
                        {
                            
                            uNormal=lowerBound+rndGen.sample01<scalar>()*(upperBound-lowerBound);

                        } while ( (sCosTheta-uNormal)/(sCosTheta-uNormalMax)*exp(sqr(uNormalMax)-sqr(uNormal)) < rndGen.sample01<scalar>());
                    }
                    else
                    {
                        uNormal = -sqrt(-log(rndGen.sample01<scalar>()));
                    }

                    U = rndGen.GaussNormal<scalar>()/sqrt(2.0)*t1
                        + rndGen.GaussNormal<scalar>()/sqrt(2.0)*t2
                        - uNormal*n;

                    gamma=1.0+(2.0/mostProbableSpeed*(faceHeatFlux.x()*U.x()+faceHeatFlux.y()*U.y()+faceHeatFlux.z()*U.z())*(2.0/5.0*(sqr(U.x())+sqr(U.y())+sqr(U.z()))-1.0)
                        -2.0*(faceStress.xy()*U.x()*U.y()+faceStress.xz()*U.x()*U.z()+faceStress.yz()*U.y()*U.z())
                        -faceStress.xx()*(sqr(U.x())-sqr(U.z()))-faceStress.yy()*(sqr(U.y())-sqr(U.z())))/pressure;

                } while(amplitudeParameter*rndGen.sample01<scalar>()>gamma);

                U=mostProbableSpeed*U+faceVelocity;
                /////////////////////////////////////////////////////////////////////////

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceRotationalTemperature,
                    cloud_.constProps(typeId).rotationalDoF()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    faceVibrationalTemperature,
                    cloud_.constProps(typeId).vibrationalDoF(),
                    typeId
                );

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceElectronicTemperature,
                    cloud_.constProps(typeId).degeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList(),
                    typeId
                );


                label newParcel = 1;

                scalar CWF = cloud_.cellWF(cellI);
                scalar RWF = cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);
                
                // Apply tracking correction towards cell centre
                p += VSMALL*(mesh_.cellCentres()[cellI] - p);

                cloud_.addNewParcel
                (
                    p,
                    U,
                    CWF,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    typeId,
                    newParcel,
                    vibLevel
                );

            }
        }
    }
}

void Foam::uniGasGeneralBoundary::insertParcels
(
    const scalar& transT,
    const vectorField& velocity
)
{

    Random& rndGen = cloud_.rndGen();

    // Loop over all species
    forAll(accumulatedParcelsToInsert_, iD)
    {
        // Loop over all faces of the patch
        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const vector& faceVelocity = velocity[f];
            const scalar faceTemperature = transT;
            const label faceI = faces_[f];
            const label cellI = cells_[f];
            const vector& fC = mesh_.faceCentres()[faceI];
            const vector& sF = mesh_.faceAreas()[faces_[f]];
            scalar fA = mag(sF);

            List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
            (
                mesh_,
                faceI,
                cellI
            );

            // Cumulative triangle area fractions
            List<scalar> cTriAFracs(faceTets.size(), 0.0);

            scalar previousCummulativeSum = 0.0;

            forAll(faceTets, triI)
            {
                const tetIndices& faceTetIs = faceTets[triI];

                cTriAFracs[triI] =
                    faceTetIs.faceTri(mesh_).mag()/fA
                    + previousCummulativeSum;

                previousCummulativeSum = cTriAFracs[triI];
            }

            // Force the last area fraction value to 1.0 to avoid any
            // rounding/non - flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            // Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = sF;
            n /= -mag(n);

            //  Wall tangential unit vector. Use the direction between the
            // face centre and the first vertex in the list
            vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
            t1 /= mag(t1);

            // Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            /* -------------------------------------------------------------*/

            scalar& faceAccumulator = accumulatedParcelsToInsert_[iD][f];

            // Number of whole particles to insert
            label nParcelsToInsert = max(label(faceAccumulator), 0);
            
            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of faceAccumulator
            if ((faceAccumulator - nParcelsToInsert) > rndGen.sample01<scalar>())
            {
                ++nParcelsToInsert;
            }

            faceAccumulator -= nParcelsToInsert;

            const label typeId = typeIds_[iD];
            const scalar pMass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nParcelsToInsert; i++)
            {
                // Choose a triangle to insert on, based on their relative
                // area

                scalar triSelection = rndGen.sample01<scalar>();

                // Selected triangle
                label selectedTriI = -1;

                forAll(cTriAFracs, triI)
                {
                    selectedTriI = triI;

                    if (cTriAFracs[triI] >= triSelection)
                    {
                        break;
                    }
                }
                
                // Randomly distribute the points on the triangle.

                const tetIndices& faceTetIs = faceTets[selectedTriI];

                point pt(faceTetIs.faceTri(mesh_).randomPoint(rndGen));

                // Velocity generation
                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        faceTemperature,
                        pMass
                    )
                );

                scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                // Coefficients required for Bird eqn 12.5
                scalar uNormProbCoeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                scalar uNormProbCoeffB =
                    0.5*
                    (
                        1.0
                        + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                    );

                // Equivalent to the QA value in Bird's DSMC3.FOR
                scalar randomScaling = 3.0;

                if (sCosTheta < -3)
                {
                    randomScaling = mag(sCosTheta) + 1;
                }

                scalar P = -1;

                // Normalised candidates for the normal direction velocity
                // component
                scalar uNormal;
                scalar uNormalThermal;

                if (mag(faceVelocity & n) > VSMALL)
                {
                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                                *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.sample01<scalar>());
                }
                else
                {
                    uNormal = sqrt(-log(rndGen.sample01<scalar>()));
                }

                vector U =
                    sqrt(physicoChemical::k.value()*faceTemperature/pMass)
                   *(
                        rndGen.GaussNormal<scalar>()*t1
                        + rndGen.GaussNormal<scalar>()*t2
                    )
                  + (t1 & faceVelocity)*t1
                  + (t2 & faceVelocity)*t2
                  + mostProbableSpeed*uNormal*n;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceTemperature,
                    cloud_.constProps(typeId).rotationalDoF()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).vibrationalDoF(),
                    typeId
                );

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).degeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList(),
                    typeId
                );
                
                label newParcel = 1;

                scalar CWF = cloud_.cellWF(cellI);
                scalar RWF = cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);
                
                pt += (n*SMALL);

                cloud_.addNewParcel
                (
                    pt,
                    U,
                    CWF,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    typeId,
                    newParcel,
                    vibLevel
                );

            }
        }
    }
}

void Foam::uniGasGeneralBoundary::insertParcels
(
    const scalarField& transT,
    const vectorField& velocity
)
{
    Random& rndGen = cloud_.rndGen();

    // loop over all species
    forAll(accumulatedParcelsToInsert_, iD)
    {
        // loop over all faces of the patch
        forAll(accumulatedParcelsToInsert_[iD], f)
        {
            const vector& faceVelocity = velocity[f];
            const scalar faceTemperature = transT[f];
            const label faceI = faces_[f];
            const label cellI = cells_[f];
            const vector& fC = mesh_.faceCentres()[faceI];
            const vector& sF = mesh_.faceAreas()[faces_[f]];
            scalar fA = mag(sF);

            List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
            (
                mesh_,
                faceI,
                cellI
            );

            //Cumulative triangle area fractions
            List<scalar> cTriAFracs(faceTets.size(), 0.0);

            scalar previousCummulativeSum = 0.0;

            forAll(faceTets, triI)
            {
                const tetIndices& faceTetIs = faceTets[triI];

                cTriAFracs[triI] =
                    faceTetIs.faceTri(mesh_).mag()/fA
                    + previousCummulativeSum;

                previousCummulativeSum = cTriAFracs[triI];
            }


            //Force the last area fraction value to 1.0 to avoid any
            //rounding/non-flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            //Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = sF;
            n /= -mag(n);

            //Wall tangential unit vector. Use the direction between the
            //face centre and the first vertex in the list
            vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
            t1 /= mag(t1);

            //Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            label nParcelsToInsert = label(accumulatedParcelsToInsert_[iD][f]);

            if ((nParcelsToInsert - accumulatedParcelsToInsert_[iD][f]) > rndGen.sample01<scalar>())
            {
                ++nParcelsToInsert;
            }

            // Note: remainder has been set
            accumulatedParcelsToInsert_[iD][f] -= nParcelsToInsert;

            const label typeId = typeIds_[iD];
            const scalar mass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nParcelsToInsert; ++i)
            {
                // Choose a triangle to insert on, based on their relative
                // area

                scalar triSelection = rndGen.sample01<scalar>();

                // Selected triangle
                label selectedTriI = -1;

                forAll(cTriAFracs, triI)
                {
                    selectedTriI = triI;

                    if (cTriAFracs[triI] >= triSelection)
                    {
                        break;
                    }
                }

                // Randomly distribute the points on the triangle.

                const tetIndices& faceTetIs = faceTets[selectedTriI];

                point p (faceTetIs.faceTri(mesh_).randomPoint(rndGen));

                // Velocity generation
                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        faceTemperature,
                        mass
                    )
                );

                scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

                // Coefficients required for Bird eqn 12.5
                scalar uNormProbCoeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                scalar uNormProbCoeffB =
                    0.5
                   *(
                        1.0
                      + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                    );

                // Equivalent to the QA value in Bird's DSMC3.FOR
                scalar randomScaling = 3.0;

                if (sCosTheta < -3)
                {
                    randomScaling = mag(sCosTheta) + 1;
                }

                scalar P = -1;

                // Normalised candidates for the normal direction velocity
                // component
                scalar uNormal;
                scalar uNormalThermal;

                if (mag(faceVelocity & n) > VSMALL)
                {
                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                               *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.sample01<scalar>());
                }
                else
                {
                    uNormal = sqrt(-log(rndGen.sample01<scalar>()));
                }

                vector U =
                    sqrt(physicoChemical::k.value()*faceTemperature/mass)
                   *(
                        rndGen.GaussNormal<scalar>()*t1
                      + rndGen.GaussNormal<scalar>()*t2
                    )
                  + (t1 & faceVelocity)*t1
                  + (t2 & faceVelocity)*t2
                  + mostProbableSpeed*uNormal*n;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceTemperature,
                    cloud_.constProps(typeId).rotationalDoF()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).vibrationalDoF(),
                    typeId
                );

                label newParcel = 1;

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).degeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList(),
                    typeId
                );

                scalar CWF = cloud_.cellWF(cellI);
                scalar RWF = cloud_.axiRWF(cloud_.mesh().cellCentres()[cellI]);
                
                // Apply tracking correction towards cell centre
                p += VSMALL*(mesh_.cellCentres()[cellI] - p);

                cloud_.addNewParcel
                (
                    p,
                    U,
                    CWF,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    typeId,
                    newParcel,
                    vibLevel
                );

            }
        }
    }
    
}


const Foam::labelList& Foam::uniGasGeneralBoundary::controlPatch() const
{
    return faces_;
}


const Foam::labelList& Foam::uniGasGeneralBoundary::controlZone() const
{
    return cells_;
}


// ************************************************************************* //
