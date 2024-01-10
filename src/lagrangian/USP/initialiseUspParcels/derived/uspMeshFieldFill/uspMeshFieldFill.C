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

#include "uspMeshFieldFill.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(uspMeshFieldFill, 0);

addToRunTimeSelectionTable
(
    uspConfiguration,
    uspMeshFieldFill,
    dictionary
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uspMeshFieldFill::uspMeshFieldFill(uspCloud& cloud, const dictionary& dict)
:
    uspConfiguration(cloud, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::uspMeshFieldFill::setInitialConfiguration()
{

    Info<< nl << "Initialising particles" << endl;

    const word initialDistributionType
    (
        uspInitialiseDict_.get<word>("initialDistributionType")
    );

    //- read initial fields
    volScalarField initialTransT_
    (
        IOobject
        (
            "transT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volScalarField initialRotT_
    (
        IOobject
        (
            "rotT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volScalarField initialVibT_
    (
        IOobject
        (
            "vibT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volScalarField initialElecT_
    (
        IOobject
        (
            "elecT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volVectorField initialU_
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );

    volVectorField initialHeatFlux_
    (
        IOobject
        (
            "heatFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimPressure*dimVelocity, Zero)
    );

    volTensorField initialStress_
    (
        IOobject
        (
            "stress",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor(dimPressure, Zero)
    );

    // Read in the type ids
    List<word> molecules = uspInitialiseDict_.lookup("typeIdList");

    //- List of inlet densities (one entry for each species)
    List<autoPtr<volScalarField>> initialNumberDensityPtr_;

    initialNumberDensityPtr_.setSize(molecules.size());

    forAll(initialNumberDensityPtr_, i)
    {

        const word& moleculeName = molecules[i];

        initialNumberDensityPtr_[i].reset
        (
            new volScalarField
            (
                IOobject
                (
                    "numberDensity_"+moleculeName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

    }

    const auto& meshCC = cloud_.mesh().cellCentres();
    const auto& meshV = cloud_.mesh().V();

    //Compute cell weights
    if (cloud_.cellWeighted())
    {
        forAll(mesh_.cells(), celli)
        {

            label geometricDims = 0;
            forAll(mesh_.geometricD(), dim)
            {
                if (mesh_.geometricD()[dim] == 1)
                {
                    geometricDims++;
                }    
            }

            scalar totalNumberDensity = 0.0;
            forAll(molecules, i)
            {
                const volScalarField& initialNumberDensity_ = initialNumberDensityPtr_[i];
                totalNumberDensity += initialNumberDensity_[celli];
            }

            scalar RWF = cloud_.axiRWF(meshCC[celli]);
            cloud_.cellWeightFactor().primitiveFieldRef()[celli] =
                (totalNumberDensity*meshV[celli])/(cloud_.particlesPerSubcell()*pow(cloud_.subcellLevels()[celli],geometricDims)*cloud_.nParticle()*RWF);

        }
        cloud_.cellWeightFactor().correctBoundaryConditions();
    }
    else
    {
        cloud_.cellWeightFactor().primitiveFieldRef() = 1.0;
        cloud_.cellWeightFactor().correctBoundaryConditions();
    }

    forAll(mesh_.cells(), celli)
    {
        List<tetIndices> cellTets =
            polyMeshTetDecomposition::cellTetIndices(mesh_, celli);

        for (const tetIndices& cellTetIs : cellTets)
        {
            tetPointRef tet = cellTetIs.tet(mesh_);

            scalar tetVolume = tet.mag();

            // compute total Pressure for Chapman-Enskog distribution
            scalar pressure=0e0;
            forAll(molecules, i)
            {
                const volScalarField& initialNumberDensity_ = initialNumberDensityPtr_[i];
                pressure += initialNumberDensity_[celli]*physicoChemical::k.value()*initialTransT_[celli];
            }

            forAll(molecules, i)
            {
                const word& moleculeName = molecules[i];

                label typeId(cloud_.typeIdList().find(moleculeName));

                if (typeId == -1)
                {
                    FatalErrorInFunction
                        << "typeId " << moleculeName << "not defined."
                        << abort(FatalError);
                }

                const auto& cP = cloud_.constProps(typeId);

                const volScalarField& initialNumberDensity_ = initialNumberDensityPtr_[i];
                scalar numberDensity = initialNumberDensity_[celli]/cloud_.nParticle();

                scalar translationalTemperature = initialTransT_[celli];

                scalar rotationalTemperature = initialRotT_[celli];

                scalar vibrationalTemperature = initialVibT_[celli];

                scalar electronicTemperature = initialElecT_[celli];

                vector velocity = initialU_[celli];

                // Calculate the number of particles required
                scalar CWF = cloud_.cellWF(celli);
                scalar RWF = cloud_.axiRWF(meshCC[celli]);
                scalar particlesRequired = numberDensity*tetVolume/(CWF*RWF);

                // Only integer numbers of particles can be inserted
                label nParticlesToInsert = label(particlesRequired);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of particlesRequired
                if
                (
                    (particlesRequired - nParticlesToInsert)
                  > rndGen_.sample01<scalar>()
                )
                {
                    ++nParticlesToInsert;
                }

                // Compute breakdown parameter if particles are initialized
                //based on the Chapman-Enskog distribution
                scalar breakdownParameter=0.0;
    
                vector heatFlux=vector::zero;

                tensor stress=tensor::zero;

                if (initialDistributionType == "ChapmanEnskog")
                {
                    scalar mostProbableSpeed(
                        cloud_.maxwellianMostProbableSpeed
                        (
                            translationalTemperature,
                            cP.mass()
                        )
                    );

                    heatFlux = initialHeatFlux_[celli];

                    stress = initialStress_[celli];

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

                    breakdownParameter = max(maxHeatFlux,maxStress);

                }

                for (label pI = 0; pI < nParticlesToInsert; ++pI)
                {
                    point p = tet.randomPoint(rndGen_);

                    vector U=vector::zero;
                    if (initialDistributionType == "Maxwellian")
                    {

                        U = cloud_.equipartitionLinearVelocity
                        (
                            translationalTemperature,
                            cP.mass()
                        );
                        
                    }
                    else if (initialDistributionType == "ChapmanEnskog")
                    {
                        U = cloud_.chapmanEnskogLinearVelocity
                        (
                            breakdownParameter,
                            pressure,
                            translationalTemperature,
                            cP.mass(),
                            heatFlux,
                            stress
                        );

                    }
                    else
                    {
                        FatalErrorIn("Foam::uspCloud<uspParcel>::initialise")
                            << "initialDistributionType " << initialDistributionType << " not defined." << nl
                            << abort(FatalError);
                    }


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

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the maximum most probable thermal speed (Bird,
    // p222 - 223)

    List<scalar> numberDensities_;
    numberDensities_.setSize(molecules.size());

    forAll(mesh_.cells(), celli)
    {

        forAll(molecules, i)
        {
            const volScalarField& initialNumberDensity_ = initialNumberDensityPtr_[i];
            numberDensities_[i] = initialNumberDensity_[celli];
        }
        
        label mostAbundantType(findMax(numberDensities_));

        const auto& cP = cloud_.constProps(mostAbundantType);

        cloud_.sigmaTcRMax()[celli] =
        cP.sigmaT()
        *cloud_.maxwellianMostProbableSpeed
            (
                initialTransT_[celli],
                cP.mass()
            );

    }

    cloud_.sigmaTcRMax().correctBoundaryConditions();

}


// ************************************************************************* //
