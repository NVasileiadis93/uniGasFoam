# About uniGasFoam usage
Using uniGasFoam will be relatively straightforward for users who are already familiar with OpenFOAM solvers. The uniGasFoam solver is based on the well-established dsmcFoam+ solver and thus, share a lot of similarities.

# Case file structure
The uniGasFoam solver follows the typical file structure of any OpenFOAM application. The base directory (referred to as *[case]*) is contains all the directories and files required to run to run a uniGasFoam simulation. Under this directory, two additional directories named *system* and *constant* are needed.

The majority of the OpenFOAM dictionary files, which control most of the running parameters for the uniGasFoam case (such as time-step control, case initialization parameters, and boundary conditions), are contained in the *[case]/system* directory. The *[case]/constant* directory contains the mesh files in [case]/constant/polyMesh* directory, as well as, the constant gas and simulation properties.

Although some dictionaries within this structure are specific to uniGasFoam, they are named using typical nomenclature and are designed to be self-descriptive.

# Dictionaries
A list of the required and optional dictionaries along with their description are given here:

**[case]/constant/**
* **uniGasProperties:** contains the constant gas and numerical simulation properties

**[case]/system/**
* **blockMeshDict**: contains data for mesh creatin with the blockMesh utility
* **controlDict**: time-control properties
* **fvSolution** (not used): numerical solution properties
* **fvSchemes**: numerical schemes and associated parameters
* **fvOptions**  (not used): solver options
* **uniGasInitialiseDict**: particle initialization method and parameters
* **boundariesDict**: boundary conditions and associated parameters
* **fieldPropertiesDict**: desired outputs 
* **chemReactDict**: chemical reaction data
* **hybridDecompositionDict** (required for hybrid schemes): continuum-rarefied domain decomposition method and associated parameters 
* **decomposeParDict** (required for parallel simulations) : spatial decompotion for parallel cases

It is noted that all numerical parameters required are always given in SI units.

# uniGasProperties dictionary
An example of the uniGasProperties dictionary is given here:
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2306                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      uniGasProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solution
{
    active          yes;
    coupled         yes;
    cellValueSourceCorrection no;
    calcFrequency   1;
    maxTrackTime    1e9;
    sourceTerms
    {
        resetOnStartup  no;
        schemes
        {}
    }
}

// General Properties
// ~~~~~~~~~~~~~~~~~~
nEquivalentParticles            5.85E+11;
chemicalReactions				false;
chargedParticles                false;
adsorption                      false;

// Cell Weighting Properties
// ~~~~~~~~~~~~~~~~~~
cellWeightedSimulation      	true;
cellWeightedProperties
{
    minParticlesPerSubCell		20;
    particlesPerSubCell         20;
}

// Axisymmetric Properties
// ~~~~~~~~~~~~~~~~~~
axisymmetricSimulation          true;
axisymmetricProperties
{
    radialExtentOfDomain        7.5e-2;
    maxRadialWeightingFactor    1000;
}

// Dynamic Adaptation
// ~~~~~~~~~~~~~~~~~~
adaptiveSimulation		     	true;
adaptiveProperties
{
    timeStepAdaptation		    true;
	subCellAdaptation		    true;
	cellWeightAdaptation		true;
	adaptationInterval			10;
    maxTimeStepMCTRatio         0.2;
    maxCourantNumber            0.5;
	maxSubCellSizeMFPRatio		1.0;
}

// Collision Models
// ~~~~~~~~~~~~~~~~~~~~~~
collisionModel              hybrid; //binary, relaxation, hybrid
relaxationCollisionModel    USBGK;
binaryCollisionPartnerModel noTimeCounter;
binaryCollisionModel        variableHardSphere;
collisionProperties
{
    Tref                    273;
    macroInterpolation      false;
    theta                   1e-3;
}

// Chemical Reactions Model
//~~~~~~~~~~~~~~~~~~~~~~~~~
//ChemicalReactionModel          QuantumKineticModel;

// Molecular species
// ~~~~~~~~~~~~~~~~~
typeIdList  (Ar);

moleculeProperties
{
  	Ar
    {
        mass                            		66.3e-27;
        diameter                        		4.17E-10;
        rotationalDegreesOfFreedom        		0;
		vibrationalModes        				0;
        omega                           		0.81;
        alpha                                   1.0;
		characteristicVibrationalTemperature	();
        dissociationTemperature            		();
        ionisationTemperature            		0;
		charDissQuantumLevel					();
		Zref 									();
		referenceTempForZref                    ();
        charge                                  0;
        numberOfElectronicLevels                1;
        electronicEnergyList                    (0);
        degeneracyList                          (1);
    }
}

// ************************************************************************* //
```

Referring to the example above, lines 1-16 are the standard OpenFOAM dictionary header and must be included in all dictionaries. 

**solution**: defines the particle tracking algorithm and should be left as is

**nEquivalentParticles**: defines the number of real gas particles that each computational particle represents

**chemicalReactions**: boolean parameter for chemical reactions, only words for pure DSMC simulations

**chargedParticles**: boolean parameter for particle ionization, only words for pure DSMC simulations

**chargedParticles**: CHECK!!!

**cellWeightedSimulation**: boolean parameter for non-uniform cell weighting 

**minParticlesPerSubCell**: minimum allowed particles per sub-cell

**particlesPerSubCell**: desired number of particles per sub-cell

**axisymmetricSimulation**: boolean parameter for axisymmetric simulation

**radialExtentOfDomain**: maximum domain radius

**maxRadialWeightingFactor**: weighting factor at maximum domain radius

**adaptiveSimulation**: boolean parameter for adaptive schemes

**timeStepAdaptation**: boolean parameter for adaptive global time steping algorithm

**subCellAdaptation**: boolean parameter for transient adaptive sub-cells algorithm

**cellWeightAdaptation**: boolean parameter for adaptive non-uniform cell weighting

**adaptationInterval**: time interval in number of time steps when adaptive properties are recalculated

**maxTimeStepMCTRatio**: maximum time step to local mean collision time ratio (reccomended <0.2)

**maxCourantNumber**: maximum Courant number (reccomended <0.5)

**maxSubCellSizeMFPRatio**: maximum subcell to local mean free path ratio (reccomended <0.5)

**collisionModel**: collision model for pure DSMC (binary), pure SP/USP (relaxation) and hybrid SP/USP (hybrid) gas flow simulation

**relaxationCollisionModel**: SP/USP relaxation model (BGK,SBGK,ESBGK,USBGK)

**binaryCollisionPartnerModel**: DSMC binary collision partner selection (noTimeCounter, noTimeCounterSubCycled)

**binaryCollisionModel**: DSMC binary collision model (variableHardSphere, variableSoftScphere) CHECK!!!

**Tref**: reference particle diameter reference temperature

**macroInterpolation**: boolean parameter for spatial interpolation of SP/USP macroscopic properties (required definition of interpolation schemes in fvSchems dictionary)

**theta**: time-averaging coefficient for SBGK and USBGK relaxation models

**ChemicalReactionModel**: chemical reaction model ()

**typeIdList**: list of gas species

**moleculeProperties**: sub-dictionary containing the properties of all the gas species defined in typeIdList

# blockMeshDict dictionary
