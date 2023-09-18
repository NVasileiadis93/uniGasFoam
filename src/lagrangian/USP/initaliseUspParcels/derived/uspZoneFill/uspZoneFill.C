/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "uspZoneFill.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspZoneFill, 0);

addToRunTimeSelectionTable
(
    uspConfiguration,
    uspZoneFill,
    dictionary
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspZoneFill::uspZoneFill(uspCloud& cloud, const dictionary& dict)
:
    uspConfiguration(cloud, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uspZoneFill::setInitialConfiguration()
{
    Info<< nl << "Initialising particles" << endl;

    const scalar translationalTemperature
    (
        uspInitialiseDict_.get<scalar>("translationalTemperature")
    );

    const scalar rotationalTemperature
    (
        uspInitialiseDict_.get<scalar>("rotationalTemperature")
    );

    const scalar vibrationalTemperature
    (
        uspInitialiseDict_.get<scalar>("vibrationalTemperature")
    );

    const scalar electronicTemperature
    (
        uspInitialiseDict_.get<scalar>("electronicTemperature")
    );

    const vector velocity(uspInitialiseDict_.get<vector>("velocity"));

    const dictionary& numberDensitiesDict =
        uspInitialiseDict_.subDict("numberDensities");

    List<word> molecules(numberDensitiesDict.toc());

    Field<scalar>numberDensities(molecules.size());

    forAll(molecules, i)
    {
        numberDensities[i] = numberDensitiesDict.get<scalar>(molecules[i]);
    }

    numberDensities /= cloud_.nParticle();

    const cellZoneMesh& cellZones = mesh_.cellZones();
    const word zoneName(uspInitialiseDict_.get<word>("zone"));
    const label zoneId = cellZones.findZoneID(zoneName);

    if (zoneId == -1)
    {
        FatalErrorInFunction
            << "Cannot find region: " << zoneName << nl << "in: "
            << mesh_.time().system()/"uspInitialiseDict"
            << exit(FatalError);
    }

    const cellZone& zone = cellZones[zoneId];

    const auto& meshCC = cloud_.mesh().cellCentres();
    const auto& meshV = cloud_.mesh().V();

    //Compute cell weights
    if (cloud_.cellWeighted())
    {
        if (zone.size())
        {

            for (const label celli : zone)
            {

                scalar totalNumberDensity = 0.0;
                forAll(molecules, i)
                {
                    totalNumberDensity += numberDensities[i]*cloud_.nParticle();;
                }

                scalar RWF = cloud_.axiRWF(meshCC[celli]);
                cloud_.cellWeightFactor().primitiveFieldRef()[celli] =
                    (totalNumberDensity*meshV[celli])/(cloud_.particlesPerCell()*cloud_.nParticle()*RWF);

            }
            cloud_.cellWeightFactor().correctBoundaryConditions();

        }

    }
    else
    {
        cloud_.cellWeightFactor().primitiveFieldRef() = 1.0;
        cloud_.cellWeightFactor().correctBoundaryConditions();
    }

    if (zone.size())
    {
        Info << "Lattice in zone: " << zoneName << endl;

        for (const label celli : zone)
        {
            List<tetIndices> cellTets =
                polyMeshTetDecomposition::cellTetIndices(mesh_, celli);

            for (const tetIndices& cellTetIs : cellTets)
            {
                tetPointRef tet = cellTetIs.tet(mesh_);

                scalar tetVolume = tet.mag();

                forAll(molecules, i)
                {
                    const word& moleculeName = molecules[i];

                    label typeId(cloud_.typeIdList().find(moleculeName));

                    if (typeId == -1)
                    {
                        FatalErrorInFunction
                            << "typeId " << moleculeName << "not defined." << nl
                            << abort(FatalError);
                    }

                    const auto& cP = cloud_.constProps(typeId);

                    scalar numberDensity = numberDensities[i];

                    // Calculate the number of particles required
                    scalar CWF = cloud_.cellWF(celli);
                    scalar RWF = cloud_.axiRWF(meshCC[celli]);
                    scalar particlesRequired = numberDensity*tetVolume/(CWF*RWF);

                    // Only integer numbers of particles can be inserted
                    label nParticlesToInsert = label(particlesRequired);

                    // Add another particle with a probability proportional to
                    // the remainder of taking the integer part of
                    // particlesRequired
                    if
                    (
                        (particlesRequired - nParticlesToInsert)
                      > rndGen_.sample01<scalar>()
                    )
                    {
                        ++nParticlesToInsert;
                    }

                    for (label pI = 0; pI < nParticlesToInsert; ++pI)
                    {
                        point p = tet.randomPoint(rndGen_);

                        vector U = cloud_.equipartitionLinearVelocity
                        (
                            translationalTemperature,
                            cP.mass()
                        );

                        scalar ERot = cloud_.equipartitionRotationalEnergy
                        (
                            rotationalTemperature,
                            cP.rotationalDoF()
                        );

                        labelList vibLevel =
                            cloud_.equipartitionVibrationalEnergyLevel
                            (
                                vibrationalTemperature,
                                cP.vibrationalDoF(),
                                typeId
                            );

                        label ELevel = cloud_.equipartitionElectronicLevel
                        (
                            electronicTemperature,
                            cP.degeneracyList(),
                            cP.electronicEnergyList(),
                            typeId
                        );

                        U += velocity;

                        label newParcel = 0;

                        scalar CWF = cloud_.cellWF(celli);
                        scalar RWF = cloud_.axiRWF(meshCC[celli]);

                        cloud_.addNewParcel
                        (
                            p,
                            U,
                            CWF,
                            RWF,
                            ERot,
                            ELevel,
                            celli,
                            typeId,
                            newParcel,
                            vibLevel
                        );
                    }
                }
            }
        }
    }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222 - 223)

    label mostAbundantType(findMax(numberDensities));

    const auto& cP = cloud_.constProps(mostAbundantType);

    for (const label celli : zone)
    {
        cloud_.sigmaTcRMax().primitiveFieldRef()[celli] =
            cP.sigmaT()*cloud_.maxwellianMostProbableSpeed
            (
                translationalTemperature,
                cP.mass()
            );
    }

    cloud_.sigmaTcRMax().correctBoundaryConditions();
}


// ************************************************************************* //
