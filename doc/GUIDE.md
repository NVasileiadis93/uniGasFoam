# About uniGasFoam usage
Using uniGasFoam will be relatively straightforward for users who are already familiar with OpenFOAM solvers. The uniGasFoam solver is based on the well-established dsmcFoam+ solver and thus, share a lot of similarities.

# Running uniGasFoam
The tutorial cases found in the *tutorials/uniGasFoam* directory are a great starting point to run and test uniGasFoam.

To run uniGasFoam as a background process and output the case progress to a log file, execute:
* uniGasFoam > log &

Like any other OpenFOAM solver, uniGasFoam can be run in parallel by decomposing the flow domain. To decompose the case and run uniGasFoam in parallel in the background and output the case progress to a log file, execute in sequence:
* decomposePar
* uniGasFoam -n [nprocs] -parallel > log &

# Case file structure
The uniGasFoam solver follows the typical file structure of any OpenFOAM application. The base directory (referred to as *[case]*) contains all the directories and files required to run a uniGasFoam simulation. Under this directory, two additional directories named *system* and *constant* are needed.

The majority of the OpenFOAM dictionary files, which control most of the running parameters for the uniGasFoam case (such as time-step control, case initialization parameters, and boundary conditions), are contained in the *[case]/system* directory. The *[case]/constant* directory contains the mesh files in [case]/constant/polyMesh* directory, as well as, the constant gas and simulation properties.

Although some dictionaries within this structure are specific to uniGasFoam, they are named using typical nomenclature and are designed to be self-descriptive.

# Dictionaries
A list of the required and optional dictionaries along with their description are given here:

**[case]/constant/**
* **uniGasProperties:** contains the constant gas and numerical simulation properties.

**[case]/system/**
* **blockMeshDict**: contains data for creating the mesh with the blockMesh utility.
* **controlDict**: time-control properties.
* **fvSolution**: numerical solution properties, not used but must exist.
* **fvSchemes**: numerical schemes and associated parameters.
* **fvOptions**: solver options, not used but must exist.
* **uniGasInitialisationDict**: particle initialisation method and parameters.
* **boundariesDict**: boundary conditions and associated parameters.
* **fieldPropertiesDict**: desired outputs.
* **chemReactDict**: chemical reaction data.
* **hybridDecompositionDict**: continuum-rarefied domain decomposition method and associated parameters, only required for hybrid simulations.
* **decomposeParDict**: multi processor spatial decomposition, only required for parallel simulations.

It is noted that all numerical parameters required are always given in SI units.

# uniGasProperties

<details>
<summary>Click to expand/collapse details</summary>
<p>

An example of the uniGasProperties dictionary is given here:
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
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
collisionModel              hybrid;
bgkCollisionModel           unifiedStochasticParticleSBGK;
dsmcCollisionPartnerModel   noTimeCounter;
dsmcCollisionModel          variableHardSphere;
collisionProperties
{
    Tref                    273;
    macroInterpolation      false;
    theta                   1e-3;
}

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

**solution**: defines the particle tracking algorithm and should be left as is.\
**nEquivalentParticles**: defines the number of real gas particles that each computational particle represents.\
**chemicalReactions**: boolean parameter for chemical reactions, only works for pure DSMC simulations.\
**chargedParticles**: boolean parameter for particle ionization, only works for pure DSMC simulations.\
**cellWeightedSimulation**: boolean parameter for non-uniform cell weighting.\
**minParticlesPerSubCell**: minimum allowed particles per sub-cell.\
**particlesPerSubCell**: desired number of particles per sub-cell.\
**axisymmetricSimulation**: boolean parameter for axisymmetric simulation.\
**radialExtentOfDomain**: maximum domain radius.\
**maxRadialWeightingFactor**: weighting factor at maximum domain radius.\
**adaptiveSimulation**: boolean parameter for adaptive schemes.\
**timeStepAdaptation**: boolean parameter for adaptive global time stepping algorithm.\
**subCellAdaptation**: boolean parameter for transient adaptive sub-cells algorithm.\
**cellWeightAdaptation**: boolean parameter for adaptive non-uniform cell weighting.\
**adaptationInterval**: time interval in number of time steps when adaptive properties are recalculated.\
**maxTimeStepMCTRatio**: maximum time step to local mean collision time ratio (recommended < 0.2).\
**maxCourantNumber**: maximum Courant number (recommended < 0.5).\
**maxSubCellSizeMFPRatio**: maximum subcell to local mean free path ratio (recommended <0 .5).\
**collisionModel**: collision model for pure DSMC (dsmc), pure SP/USP (bgk) and hybrid SP/USP-DSMC (hybrid) gas flow simulation.\
**bgkCollisionModel**: SP/USP collision model (noBGKCollision, stochasticParticleBGK, stochasticParticleSBGK, stochasticParticleESBGK, unifiedStochasticParticleSBGK).\
**dsmcCollisionPartnerModel**: DSMC collision partner selection (noTimeCounter, noTimeCounterSubCycled).\
**dsmcCollisionModel**: DSMC collision model (noDsmcCollision, variableHardSphere, variableSoftSphere, LarsenBorgnakkeVariableSoftSphere).\
**Tref**: reference particle diameter reference temperature.\
**macroInterpolation**: boolean parameter for spatial interpolation of SP/USP macroscopic properties (required definition of interpolation schemes in fvSchems dictionary).\
**theta**: time-averaging coefficient for stochasticParticleSBGK and unifiedStochasticParticleSBGK collision models.\
**typeIdList**: list of gas species.\
**moleculeProperties**: sub-dictionary containing the properties of all the gas species defined in typeIdList.

</p>
</details>

# blockMeshDict

<details>
<summary>Click to expand/collapse details</summary>
<p>

An example of the blockMeshDict dictionary is given here:
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version         2.0;
    format          ascii;

    root            "";
    case            "";
    instance        "";
    local           "";

    class           dictionary;
    object          blockMeshDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale 0.3048;

vertices
(
    (-2.00000000000000E+00	0.00000000000000E+00	 5.00000000000000E-02)
    (-5.00000000000000E-01	0.00000000000000E+00	 5.00000000000000E-02)
    (0.00000000000000E+00	5.00000000000000E-01	 5.00000000000000E-02)
    (5.00000000000000E-01	0.00000000000000E+00	 5.00000000000000E-02)
    (2.00000000000000E+00	0.00000000000000E+00	 5.00000000000000E-02)
    (0.00000000000000E+00	2.00000000000000E+00	 5.00000000000000E-02)
    (-2.00000000000000E+00	0.00000000000000E+00	-5.00000000000000E-02)
    (-5.00000000000000E-01	0.00000000000000E+00	-5.00000000000000E-02)
    (0.00000000000000E+00	5.00000000000000E-01	-5.00000000000000E-02)
    (5.00000000000000E-01	0.00000000000000E+00	-5.00000000000000E-02)
    (2.00000000000000E+00	0.00000000000000E+00	-5.00000000000000E-02)
    (0.00000000000000E+00	2.00000000000000E+00	-5.00000000000000E-02)

);

edges
(
    arc 1	2	(-3.53553390593274E-01	3.53553390593274E-01	 5.00000000000000E-02)
    arc 2	3	(3.53553390593274E-01	3.53553390593274E-01	 5.00000000000000E-02)
    arc 0	5	(-1.41421356237310E+00	1.41421356237309E+00	 5.00000000000000E-02)
    arc 5	4	(1.41421356237310E+00	1.41421356237309E+00	 5.00000000000000E-02)
    arc 7	8	(-3.53553390593274E-01	3.53553390593274E-01	-5.00000000000000E-02)
    arc 8	9	(3.53553390593274E-01	3.53553390593274E-01	-5.00000000000000E-02)
    arc 6	11	(-1.41421356237310E+00	1.41421356237309E+00	-5.00000000000000E-02)
    arc 11	10	(1.41421356237310E+00	1.41421356237309E+00	-5.00000000000000E-02)
);

blocks
(
	hex (6 7 8 11 0 1 2 5) uniGasZone (100 100 1) simpleGrading (0.2 1 1)
	hex (2 5 4 3 8 11 10 9) uniGasZone (100 100 1) simpleGrading (5 1 1)
);

boundary
(

    inlet
    {
        type patch;
        faces (
			  (0 5 11 6)
              );
    }
    outlet
    {
        type patch;
        faces (
			  (4 5 11 10)
              );
    }
    cylinder
    {
        type wall;
        faces (
			  (1 2 8 7)
			  (2 3 9 8)
              );
    }
    axis
    {
		type symmetry;
        faces (
			  (0 1 7 6)
			  (3 4 10 9)
              );
    }
    empty
    {
        type empty;
        faces (
			  (0 5 2 1)
			  (2 5 4 3)
			  (6 11 8 7)
			  (8 9 10 11)
              );    
    }

);

mergePatchPairs
(
);

// ************************************************************************* //
```

**scale**: scaling factor for the vertex coordinates.\
**vertices**: definition of mesh vertices.\
**edges**: definition of mesh arc or spline edges.\
**blocks**: definition of mesh hexahedral blocks.\
**boundary**: definition of mesh boundary patches and types.\
**mergePatchPairs**: list of patches to be merged.
	
More detailed information about each entry and use can be found in https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.3-mesh-generation-with-the-blockmesh-utility.
</p>
</details>

# controlDict

<details>
<summary>Click to expand/collapse details</summary>
<p>

An example of the controlDict dictionary is given here:

```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     uniGasFoam;

nTerminalOutputs    10;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         4e-3;

deltaT          1e-7;

writeControl    adjustable;

writeInterval   1e-5;

purgeWrite      0;

writeFormat     ascii;

writePrecision  10;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable yes;

adjustTimeStep  no;

// ************************************************************************* //
```

**application**: name of solver.\
**nTerminalOutputs**: time interval in number of time steps when simulation info is printed on screen.\
**startFrom**: controller for simulation start time (firstTime,startTime,latestTime).\
**startTime**: simulation start time when startFrom startTime is selected.\
**stopAt**: controller for simulation end time (endTime,writeNow,noWriteNow,nextWrite).\
**endTime**: simulation end time when stopAt endTime is selected.\
**deltaT**: computational time step.\
**writeControl**: controller for timing of writing output to file (none, timeStep, runTime, adjustable, adjustableRunTime, clockTime, cpuTime).\
**writeInterval**: time interval for writing output to file used in conjunction with writeControl.\
**purgeWrite**: integer representing a limit on the number of time directories that are stored by overwriting time directories on a cyclic basis, set to 0 to disable the time directory limit.\
**writeFormat**: format of output files (ascii,binary).\
**writePrecision**: precision of output files.\
**writeCompression**: compression of output files (uncompressed, compressed).\
**timeFormat**: format of time directories naming (fixed, scientific, general).\
**timePrecision**: precision of time directories naming used in conjuction with timeFormat.\
**runTimeModifiable**: boolean parameter for whether dictionaries are re-read at the beginning of each time step.\
**adjustTimeStep**: boolean parameter for adjusting computational time step, always set to off since time-step adaptation is defined by uniGasProperties dictionary.

More information about each entry can be found in https://www.openfoam.com/documentation/user-guide/6-solving/6.1-time-and-data-inputoutput-control.

</p>
</details>

# fvSolution

<details>
<summary>Click to expand/collapse details</summary>
<p>

The fvSolution dictionary is not used by uniGasFoam so it should be left empty as shown below:
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
}

// ************************************************************************* //
```

</p>
</details>

# fvSchemes

<details>
<summary>Click to expand/collapse details</summary>
<p>

An example of the fvSchemes dictionary is given here:
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default         none;
}

gradSchemes
{
    default         none;
}

divSchemes
{
    default         none;
}

laplacianSchemes
{
    default         none;
}

interpolationSchemes
{
    default         	linear;
    Prandtl         	cellPoint;
    collFreq        	cellPoint;
    p              	    cellPoint;
    translationalT  	cellPoint;
    UMean               cellPoint;
    heatFluxVector	    cellPoint;
    shearStressTensor   cellPoint;
}

snGradSchemes
{
    default         none;
}

fluxRequired
{
    default         no;
}

// ************************************************************************* //
```

From the fvSchemes dictionary only the interpolationSchemes sub-dictionary is used by uniGasFoam. The default linear scheme is used in Laplacian smoothing implemented to reduce the statistical noise of macroscopic quantities, such as gas density, temperature and velocity used in adaptive schemes and domain decomposition modules. The cellPoint scheme is used to reconstruct macroscopic fields, such as gas Prandtl number, collision frequency, pressure, temperature etc to increase the spatial accuracy of SP and USP schemes. The macroscopic field reconstruction is only implemented if macroInterpolation is enabled in the uniGasProperties dictionary.

</p>
</details>

# fvOptions

<details>
<summary>Click to expand/collapse details</summary>
<p>

The fvOptions dictionary is not used by uniGasFoam so it should be left empty as shown below:
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvOptions;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
```

</p>
</details>

# uniGasInitialisationDict

<details>
<summary>Click to expand/collapse details</summary>
<p>

An example of the uniGasInitialisationDict dictionary is given here:
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      uniGasInitialisationDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

configurations
(
    configuration
    {
        type                        uniGasMeshFill;
        numberDensities             {Ar 3.416E+22;};
		translationalTemperature    300;
		rotationalTemperature       0;
		vibrationalTemperature      0;
        electronicTemperature       0;
		velocity		            (484 0 0); 
    }
);

// ************************************************************************* //
```

**configurations**: list of initialisation configurations.\
**type**: type of initilisation (uniGasMeshFill, uniGasZoneFill, uniGasMeshFieldFill).

## List of initialisation configurations

### uniGasMeshFill initialisation
Uniform initialisation of the entire flow domain. An example of the uniGasMeshFill configuration is given here:
```
configuration
{
    type                        uniGasMeshFill;
    numberDensities             {Ar 4.24700E+20;};
    translationalTemperature    200;    
    rotationalTemperature		200;	
    vibrationalTemperature		200;
    electronicTemperature       200;
    velocity 			        (2634.7 0 0);
}
```

**numberDensities**: list of initial number densities for each gas species.\
**translationalTemperature**: initial translational temperature.\
**rotationalTemperature**: initial rotational temperature.\
**vibrationalTemperature**: initial vibrational temperature.\
**electronicTemperature**: initial electronic temperature.\
**velocity**: initial velocity vector.

### uniGasMeshZoneFill initialisation
Uniform initialisation of a specific zone in the flow domain. The different initialisation zones must be defined for example by using the topoSet utility https://www.openfoam.com/documentation/guides/latest/doc/guide-meshing-topoSet.html. An example of the uniGasMeshZoneFill configuration is given here:
```
configuration
{
    type                        uniGasZoneFill;
    zone                        uniGasZone;	
    numberDensities             {Ar 4.24700E+20;};
    translationalTemperature    200;    
    rotationalTemperature		200;	
    vibrationalTemperature		200;
    electronicTemperature       200;
    velocity 			        (2634.7 0 0);
}
```

**zone**: zone name to be initialised.\
**numberDensities**: list of initial number densities for each gas species.\
**translationalTemperature**: initial translational temperature.\
**rotationalTemperature**: initial rotational temperature.\
**vibrationalTemperature**: initial vibrational temperature.\
**electronicTemperature**: initial electronic temperature.\
**velocity**: initial velocity vector.

### uniGasMeshFieldFill initialisation
Non-uniform initialisation of the entire flow domain. An example of the uniGasMeshZoneFill configuration is given here:

```
configurations
(
    configuration
    {
        type            uniGasMeshFieldFill;
        typeIdList      (Ar);
    }
);
```

**typeIdList** list of gas species initialised by the configuration.

The initial gas properties are given in separate dictionaries located in the initial time dictionary. The required initialisation dictionaries for this example are:

* **numberDensity_Ar**: initial number density of Argon (Ar) for each mesh cell.
* **transT**: initial translational temperature for each mesh cell.
* **rotT**: initial rotational temperature for each mesh cell.
* **vibT**: initial vibrational temperature for each mesh cell.
* **elecT**: initial electronic temperature for each mesh cell.
* **U**: initial velocity vector for each mesh cell.

</p>
</details>

# boundariesDict

<details>
<summary>Click to expand/collapse details</summary>
<p>

An example of the boundariesDict dictionary is given here:

```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      boundariesDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

uniGasPatchBoundaries
(

    boundary
    {
        patchBoundaryProperties
        {
            patch       cylinder;
        }

        boundaryModel   uniGasDiffuseWallPatch;

        uniGasDiffuseWallPatchProperties
        {
            velocity 	    (0 0 0);
            temperature     500;
        }
    }

	boundary
    {
        patchBoundaryProperties
        {
            patch       inlet;
        }

        boundaryModel   uniGasDeletionPatch;

        uniGasDeletionPatchProperties
        {
          allSpecies        yes;
        }
    }

	boundary
    {
        patchBoundaryProperties
        {
            patch       outlet;
        }

        boundaryModel   uniGasDeletionPatch;

        uniGasDeletionPatchProperties
        {
          allSpecies        yes;
        }
    }

);

uniGasGeneralBoundaries
(

    boundary
    {
        generalBoundaryProperties
        {
            patch       inlet;
        }

        boundaryModel   uniGasFreeStreamInflowPatch;

        uniGasFreeStreamInflowPatchProperties
        {
            translationalTemperature    200;
			rotationalTemperature		200;
			vibrationalTemperature		200;
            electronicTemperature       200;
            velocity 			        (2634.7 0 0);
			typeIds				        (Ar);
			numberDensities{Ar	        4.24700E+20;}
        }
    }

);

uniGasCyclicBoundaries
(
);

// ************************************************************************* //
```

**uniGasPatchBoundaries**: list of patch boundary definitions.\
**uniGasGeneralBoundaries**: list of general boundary definitions.\
**uniGasCyclicBoundaries**: list of cyclic boundary definitions.\
**generalBoundaryProperties**: general boundary condition properties.\
**patch**: boundary patch name definition.\
**boundaryModel**: definition of boundary condition model name.\
**[boundaryModel]Properties** specific properties required by each boundary condition.

## List of uniGasPatchBoundaries

### uniGasDiffuseWallPatch boundary condition
Fully diffuse boundary condition with uniform wall temperature and velocity. An example of the uniGasDiffuseWallPatch boundary condition definition is given here:
```
boundary
{
    patchBoundaryProperties
    {
        patch           wall;
    }

    boundaryModel   uniGasDiffuseWallPatch;

    uniGasDiffuseWallPatchProperties
    {
        temperature     300;
        velocity        (0 0 0);
    }
}
```

**temperature**: uniform wall temperature.\
**velocity**: uniform wall velocity vector.

### uniGasDiffuseWallFieldPatch boundary condition
Fully diffuse boundary condition with non-uniform wall temperature and velocity. An example of the uniGasDiffuseWallFieldPatch boundary condition definition is given here:
```
boundary
{
    patchBoundaryProperties
    {
        patch           wall;
    }

    boundaryModel   uniGasDiffuseWallFieldPatch;

    uniGasDiffuseWallFieldPatchProperties
    {
    }
}
```
The wall properties are given in separate dictionaries located in the initial time dictionary. The required boundary dictionaries for the uniGasDiffuseWallFieldPatch boundary condition are:

* **boundaryT**: boundary temperature for each patch face.
* **boundaryU**: boundary velocity vector for each patch face.

### uniGasSpecularWallPatch boundary condition
Fully specular boundary condition. An example of the uniGasSpecularWallPatch boundary condition definition is given here:
```
boundary
{
    patchBoundaryProperties
    {
        patchName                           wall;
    }

    boundaryModel   uniGasSpecularWallPatch;

    uniGasSpecularWallPatchProperties
    {
    }
}
```

### uniGasMixedDiffuseSpecularWallPatch boundary condition
Mixed diffuse-specular boundary condition with uniform wall temperature and velocity. An example of the uniGasMixedDiffuseSpecularWallPatch boundary condition definition is given here:
```
boundary
{
    patchBoundaryProperties
    {
        patch               wall;
    }

    boundaryModel   uniGasMixedDiffuseSpecularWallPatch;

    uniGasMixedDiffuseSpecularWallPatchProperties
    {
        diffuseFraction     0.9;
        temperature         300;
        velocity            (0 0 0);
    }
}
```

**diffuseFraction**: wall accommodation coefficient.\
**temperature**: uniform wall temperature.\
**velocity**: uniform wall velocity vector.

### uniGasMixedDiffuseSpecularWallFieldPatch boundary condition
Mixed diffuse-specular boundary condition with non-uniform wall temperature and velocity. An example of the uniGasMixedDiffuseSpecularWallFieldPatch boundary condition definition is given here:
```
boundary
{
    patchBoundaryProperties
    {
        patch               wall;
    }

    boundaryModel   uniGasMixedDiffuseSpecularWallFieldPatch;

    uniGasMixedDiffuseSpecularWallFieldPatchProperties
    {
        diffuseFraction     0.9;
    }
}
```

**diffuseFraction**: wall accommodation coefficient.

The wall properties are given in separate dictionaries located in the initial time dictionary. The required boundary dictionaries for the uniGasDiffuseWallFieldPatch boundary condition are:

* **boundaryT**: boundary temperature for each patch face.
* **boundaryU**: boundary velocity vector for each patch face.

### uniGasCLLWallPatch boundary condition
Cercignani-Lampis-Lord boundary condition with uniform wall temperature and velocity. An example of the uniGasCLLWallPatch boundary condition definition is given here:
```
boundary
{
    patchBoundaryProperties
    {
        patch       wall;
    }

    boundaryModel   uniGasCLLWallPatch;

    uniGasCLLWallPatchProperties
    {
        normalAccommCoeff       0.9;
        tangentialAccommCoeff   0.9;
        rotEnergyAccommCoeff    0.9;
        temperature             300;
        velocity                (0 0 0);
    }
}
```

**normalAccommCoeff**: wall normal energy accommodation coefficient.\
**tangentialAccommCoeff**: wall tangential momentum accommodation coefficient.\
**rotEnergyAccommCoeff**: wall rotationalEnergy accommodation coefficient.\
**temperature**: uniform wall temperature.\
**velocity**: uniform wall velocity vector.

### uniGasCLLWallFieldPatch boundary condition
Cercignani-Lampis-Lord boundary condition with non-uniform wall temperature and velocity. An example of the uniGasCLLWallFieldPatch boundary condition definition is given here:
```
boundary
{
    patchBoundaryProperties
    {
        patch       wall;
    }

    boundaryModel   uniGasCLLWallFieldPatch;

    uniGasCLLWallFieldPatchProperties
    {
        normalAccommCoeff       0.9;
        tangentialAccommCoeff   0.9;
        rotEnergyAccommCoeff    0.9;
    }
}
```

**normalAccommCoeff**: wall normal energy accommodation coefficient.\
**tangentialAccommCoeff**: wall tangential momentum accommodation coefficient.\
**rotEnergyAccommCoeff**: wall rotationalEnergy accommodation coefficient.

The wall properties are given in separate dictionaries located in the initial time dictionary. The required boundary dictionaries for the uniGasDiffuseWallFieldPatch boundary condition are:

* **boundaryT**: boundary temperature for each patch face.
* **boundaryU**: boundary velocity vector for each patch face.

### uniGasDeletionPatch boundary condition
Boundary condition for deleting particles exiting the domain. Should be specified for all open patches. An example of the uniGasDeletionPatch boundary condition definition is given here:
```
boundary
{
    patchBoundaryProperties
    {
        patch       outlet;
    }

    boundaryModel   uniGasDeletionPatch;

    uniGasDeletionPatchProperties
    {
    }
}
```

## List of uniGasGeneralBoundaries

### uniGasFreeStreamInflowPatch boundary condition
Maxwellian free-stream inlet boundary condition with uniform inlet conditions. An example of the uniGasFreeStreamInflowPatch boundary condition definition is given here:

```
boundary
    {
    generalBoundaryProperties
    {
        patch       flow;
    }
    boundaryModel   uniGasFreeStreamInflowPatch;
    uniGasFreeStreamInflowPatchProperties
    {
        typeIds                     (Ar);
        numberDensities             {Ar    1.0e21;}
        translationalTemperature    300;
        rotationalTemperature       300;
        vibrationalTemperature      300;
        electronicTemperature       300;
        velocity                    (100 0 0);
    }
}
```

**typeIds**: list of incoming gas species.\
**numberDensities**: list of uniform inlet number densities.\
**translationalTemperature**: uniform inlet translational temperature.\
**rotationalTemperature**: uniform inlet rotational temperature.\
**vibrationalTemperature**: uniform inlet vibrational temperature.\
**electronicTemperature**: uniform inlet electronic temperature.\
**velocity**: uniform inlet velocity vector

### uniGasFreeStreamInflowFieldPatch boundary condition
Maxwellian free-stream inlet boundary condition with non-uniform inlet conditions. An example of the uniGasFreeStreamInflowFieldPatch boundary condition definition is given here:

```
boundary
{
    generalBoundaryProperties
    {
        patch       flow;
    }
    boundaryModel   uniGasFreeStreamInflowFieldPatch;

    uniGasFreeStreamInflowFieldPatchProperties
    {
        typeIds     (Ar);
    }
}
```

**typeIds**: list of incoming gas species

The inlet properties are given in separate dictionaries located in the initial time dictionary. The required boundary dictionaries for the uniGasDiffuseWallFieldPatch boundary condition in this example are:

* **numberDensity_Ar**: boundary number density of Argon for each patch face.
* **boundaryTransT**: boundary translational temperature for each patch face.
* **boundaryRotT**: boundary rotational temperature for each patch face.
* **boundaryVibT**: boundary vibrational temperature for each patch face.
* **boundaryElecT**: boundary electronic temperature for each patch face.
* **boundaryU**: boundary velocity vector for each patch face.

### uniGasChapmanEnskogFreeStreamInflowPatch boundary condition
Chapman-Enskog free-stream inlet boundary condition with uniform inlet conditions. An example of the uniGasChapmanEnskogFreeStreamInflowPatch boundary condition definition is given here:

```
boundary
{
    generalBoundaryProperties
    {
        patch       flow;
    }

    boundaryModel   uniGasChapmanEnskogFreeStreamInflowPatch;

    uniGasChapmanEnskogFreeStreamInflowPatchProperties
    {
        typeIds                     (Ar);
        numberDensities             {Ar     1.0e21;}
        translationalTemperature    300;
        rotationalTemperature       300;
        vibrationalTemperature      300;
        electronicTemperature       300;
        velocity                    (100 0 0);
        heatFlux                    (0 0 0);
        stress                      (0.2 0 0.1
                                    0 0.1 0.3
                                    0.1 0.3 -0.3);
    }
}
```

**typeIds**: list of incoming gas species.\
**numberDensities**: list of uniform inlet number densities.\
**translationalTemperature**: uniform inlet translational temperature.\
**rotationalTemperature**: uniform inlet rotational temperature.\
**vibrationalTemperature**: uniform inlet vibrational temperature.\
**electronicTemperature**: uniform inlet electronic temperature
**velocity**: uniform inlet velocity vector.\
**heatFlux**: uniform heat flux  vector.\
**stress**: uniform stress tensor.

### uniGasChapmanEnskogFreeStreamInflowFieldPatch boundary condition
Chapman-Enskog free-stream inlet boundary condition with non-uniform inlet conditions. An example of the uniGasChapmanEnskogFreeStreamInflowFieldPatch boundary condition definition is given here:

```
boundary
{
    generalBoundaryProperties
    {
        patch       flow;
    }

    boundaryModel   uniGasChapmanEnskogFreeStreamInflowFieldPatch;

    uniGasChapmanEnskogFreeStreamInflowFieldPatchProperties
    {
        typeIds     (Ar);
    }
}
```

**typeIds**: list of incoming gas species

The inlet properties are given in separate dictionaries located in the initial time dictionary. The required boundary dictionaries for the uniGasChapmanEnskogFreeStreamInflowFieldPatch boundary condition in this example are:

* **numberDensity_Ar**: boundary number density of Argon for each patch face.
* **boundaryTransT**: boundary translational temperature for each patch face.
* **boundaryRotT**: boundary rotational temperature for each patch face.
* **boundaryVibT**: boundary vibrational temperature for each patch face.
* **boundaryElecT**: boundary electronic temperature for each patch face.
* **boundaryU**: boundary velocity vector for each patch face.
* **boundaryHeatFlux**: boundary heat flux vector for each patch face.
* **boundaryStress**: boundary stress tensor for each patch face.

### uniGasLiouFangPressureInletPatch boundary condition
Low-speed inlet pressure boundary condition based on the work of W. W. Liou and Y. C. Fang. An example of the uniGasLiouFangPressureInletPatch boundary condition definition is given here:

```
boundary
{
    generalBoundaryProperties
    {
        patch      inlet;
    }

    boundaryModel   uniGasLiouFangPressureInletPatch;

    uniGasLiouFangPressureInletPatchProperties
    {
        typeIds             (Ar);
        moleFractions       {Ar     1.0;}
        inletPressure       252000;
        inletTemperature    300;
        theta               0.01;
    }
}
```

**typeIds**: list of incoming gas species.\
**moleFractions**: list of mole fraction for each incoming gas specie.\
**inletPressure**: uniform inlet pressure.\
**inletTemperature**: uniform inlet temperature.\
**theta**: time-averaging coefficient.

### uniGasWangPressureInletPatch boundary condition
Low-speed inlet pressure boundary condition based on the work of M. Wang and Z. Li. An example of the uniGasWangPressureInletPatch boundary condition definition is given here:

```
boundary
{
    generalBoundaryProperties
    {
        patch       inlet;
    }

    boundaryModel   uniGasWangPressureInletPatch;

    uniGasWangPressureInletPatchProperties
    {
        typeIds             (Ar);
        moleFractions       {Ar      1.0;}
        inletPressure       252000;
        inletTemperature    300;
    }
}
```

**typeIds**: list of incoming gas species.\
**moleFractions**: list of mole fraction for each incoming gas specie.\
**inletPressure**: uniform inlet pressure.\
**inletTemperature**: uniform inlet temperature.

### uniGasLiouFangPressureOutletPatch boundary condition
Low-speed outlet pressure boundary condition based on the work of W. W. Liou and Y. C. Fang. An example of the uniGasLiouFangPressureOutletPatch boundary condition definition is given here:

```
boundary
{
    generalBoundaryProperties
    {
        patch       outlet;
    }

    boundaryModel   uniGasLiouFangPressureOutletPatch;

    uniGasLiouFangPressureOutletPatchProperties
    {
        typeIds             (Ar);
        moleFractions       {Ar     1.0;}
        outletPressure      100000;
    }
}
```

**typeIds**: list of incoming gas species.\
**moleFractions**: list of mole fraction for each incoming gas specie.\
**outletPressure**: uniform outlet pressure.

### uniGasMassFlowRateInletPatch boundary condition
Low speed mass flow rate inlet boundary condition based on the work of M. Lei et al. An example of the uniGasMassFlowRateInletPatch boundary condition definition is given here:

```
boundary
{
    generalBoundaryProperties
    {
        patch       inlet;
    }

    boundaryModel   uniGasMassFlowRateInletPatch;

    uniGasMassFlowRateInletPatchProperties
    {
        typeIds                     (Ar);
        moleFractions               {Ar     1.0;}
        inletTemperature            300;
        massFlowRate                1e-11;
        initialVelocity             (1 0 0);
        theta                       0.01;
    }
}
```

**typeIds**: list of incoming gas species.\
**moleFractions**: list of mole fraction for each incoming gas specie.\
**inletTemperature**: uniform inlet temperature.\
**massFlowRate**: inlet mass flow rate.\
**initialVelocity**: initial guess of inlet velocity.\
**theta**: time-averaging coefficient.

## List of uniGasCyclicBoundaries

### uniGasReflectiveParticleMembranePatch boundary condition
Cyclic reflective membrane boundary condition based on the work of Li et al. An example of the uniGasReflectiveParticleMembranePatch boundary condition definition is given here:

```
boundary
{
    generalBoundaryProperties
    {
        patch       membrane;
    }

    boundaryModel   uniGasReflectiveParticleMembranePatch;

    uniGasReflectiveParticleMembranePatchProperties
    {
        reflectionProbability       0.5;
        temperature                 300;
        velocity                    (0 0 0);
    }
}
```

**temperature**: membrane reflection probability.\
**temperature**: uniform membrane temperature.\
**velocity**: uniform membrane velocity vector.

## List of standard boundary conditions
In addition to the uniGasFoam specific boundary conditions defined in the *[case]/system/boundariesDict* dictionary several standard OpenFOAM boundary conditions can be used. The standard boundary conditions are implemented during the mesh creation process through the *[case]/system/blockMeshDict* dictionary when the blockMesh utility is used or by manually changing the *[case]/constant/polyMesh/boundary* file if third party meshing programs are used.

The standard boundary conditions that are currently available in uniGasFoam are:

* **cyclic**: cyclic boundary condition between a pair of boundaries.
* **symmetry**:  symmetry boundary condition for non-planar patches. Same as pure specular reflection.
* **symmetryPlane**:  symmetry boundary condition for planar patches. Same as pure specular reflection.
* **empty**:  empty boundary condition for reduced dimensions cases, i.e. 1-D and 2-D geometries. This condition is applied to patches whose normal is aligned to geometric directions that do not constitute solution directions.

</p>
</details>

# fieldPropertiesDict

<details>
<summary>Click to expand/collapse details</summary>
<p>

An example of the fieldPropertiesDict dictionary is given here:

```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
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

uniGasFields
(
     field
     {
         fieldModel          	        uniGasVolFields;
 
         timeProperties
         {
         	sampleInterval              1;
 	    	resetAtOutput		        on;
		    resetAtOutputUntilTime      2e-6;
         }
         
         uniGasVolFieldsProperties
         {
            field                       Ar;
         	typeIds                     (Ar);
         	measureMeanFreePath         true;
         	Tref                        273;
         	measureErrors               true;
		    averagingAcrossManyRuns     true;	
         }
     }

     field
     {
         fieldModel                     uniGasMassFluxSurface;

         timeProperties
         {
            sampleInterval              1;
            resetAtOutput               on;
            resetAtOutputUntilTime      2e-6;
         }

         uniGasMassFluxSurfaceProperties
         {
            field                       Ar;
            typeIds                     (Ar);
            faceZone                    nozzleOutlet;
            fluxDirection               (1 0 0);
            averagingAcrossManyRuns     true;
	    }
     }

     field
     {
         fieldModel          	        uniGasForceSurface;
 
         timeProperties
         {
            sampleInterval              1;
         	resetAtOutput               on;
         	resetAtOutputUntilTime  	2e-6;
         }
         
         uniGasForceSurfaceProperties
         {
            field                       Ar;
         	typeIds                     (Ar);
            patch                       targetSurface;
            averagingAcrossManyRuns     true;
         }
     } 
);

// ************************************************************************* //
```

**fieldModel**: name of field model (uniGasVolFields,uspMassFluxSurface,uniGasForceSurface).\
**timeProperties**: sub-dictionary of time properties \
**sampleInterval**: time-interval between successive samples.\
**resetAtOutput**: boolean for resetting samples every writeInterval.\
**resetAtOutputUntilTime**: steady-state time to stop sample resetting.\
**[fieldModel]Properties**: sub-dictionary of field model specific properties.

## List of fieldModels

### uniGasVolFields
Computation of macroscopic fields  at each cell and boundary face. An example of the uniGasVolFields field definition is given here:
```
field
{
    fieldModel          	        uniGasVolFields;

    timeProperties
    {
    	sampleInterval              1;
   	    resetAtOutput		        on;
	    resetAtOutputUntilTime      2e-6;
    }
    
    uniGasVolFieldsProperties
    {
        field                       Ar;
    	typeIds                     (Ar);
    	measureMeanFreePath         true;
    	Tref                        273;
    	measureErrors               true;
	    averagingAcrossManyRuns     true;	
    }
}
```

**field**: name of field.\
**typeIds**: gas species included in the macroscopic field computations.\
**measureMeanFreePath**: boolean for measuring mean free path related fields.\
**Tref**: reference temperature for mean free path calculation.\
**measureErrors**: boolean for measuring macroscopic field errors based on the work of N. Hadjiconstantinou et al.\
**averagingAcrossManyRuns**: boolean for storing average macroscopic fields when simulation is stopped. Needed for resuming the simulation without resetting macroscopic fields.

### uniGasMassFluxSurface
Computation of particle, mass and momentum flow across a surface. An example of the uniGasMassFluxSurface field definition is given here:
```
field
{
    fieldModel                     uniGasMassFluxSurface;

    timeProperties
    {
        sampleInterval              1;
        resetAtOutput               on;
        resetAtOutputUntilTime      2e-6;
    }

    uniGasMassFluxSurfaceProperties
    {
        field                       Ar;
        typeIds                     (Ar);
        faceZone                    nozzleOutlet;
        fluxDirection               (1 0 0);
        averagingAcrossManyRuns     true;
   }
}
```

**field**: name of field.\
**typeIds**: gas species included in the particle, mass and momentum flow computation.\
**faceZone**: name of face zone for flux calculations. Must be internal (not boundary) and defined during the mesh creation process for example by using the topoSet utility .\
**fluxDirection**: flux direction vector.\
**averagingAcrossManyRuns**: boolean for storing average macroscopic fields when simulation is stopped. Needed for resuming the simulation without resetting macroscopic fields.

### uniGasForceSurface
Computation of force acting on a solid boundary. An example of the uniGasForceSurface field definition is given here:
```
field
{
    fieldModel          	        uniGasForceSurface;
 
    timeProperties
    {
        sampleInterval              1;
    	resetAtOutput               on;
    	resetAtOutputUntilTime  	2e-6;
    }
    
    uniGasForceSurfaceProperties
    {
        field                       Ar;
        typeIds                     (Ar);
        patch                       targetSurface;
        averagingAcrossManyRuns     true;
    }
} 
```

**field**: name of field.\
**typeIds**: gas species included in the force computations.\
**patch**: name of patch for force calculation.\
**averagingAcrossManyRuns**: boolean for storing average macroscopic fields when simulation is stopped. Needed for resuming the simulation without resetting macroscopic fields.

</p>
</details>

# hybridDecompositionDict

<details>
<summary>Click to expand/collapse details</summary>
<p>

An example of the hybridDecompositionDict dictionary is given here:
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      hybridDecompositionProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

decompositionModel          localKnudsen;

timeProperties
{
    decompositionInterval           50;
    resetAtDecomposition            on;
    resetAtDecompositionUntilTime   1e-3;
}

localKnudsenProperties
{
    breakdownMax                    0.05;
    theta                           0.8;
    smoothingPasses	                25;
}
```

**decompositionModel**: name of continuum breakdown criterion.\
**timeProperties**: sub-dictionary of time properties.\
**decompositionInterval**: time-interval between successive decompositions.\
**resetAtDecomposition**: boolean for resetting decomposition quantities every writeInterval.\
**resetAtDecompositionUntilTime**: steady-state time to stop decomposition quantities resetting.\
**localKnudsenProperties**: sub-dictionary of continuum breakdown model properties.\
**breakdownMax**: continuum breakdown maximum threshold value.\
**theta**: time-averaging coefficient.\
**smoothingPasses**: number of passes for smoothing domain decomposition.

</p>
</details>

# decomposeParDict

<details>
<summary>Click to expand/collapse details</summary>
<p>

An example of the decomposeParDict dictionary is given here:
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2406                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains  20;

method              scotch;

// ************************************************************************* //
```

uniGasFoam uses the standard parallel domain decomposition implemented in OpenFOAM. More detailed information about each entry and its use can be found in https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.3-mesh-generation-with-the-blockmesh-utility.

</p>
</details>
