# Prerequisites
An OpenFOAM installation of v2212 or later is required. The latest OpenFOAM version along with detailed installation instructions can be found at:

https://develop.openfoam.com/Development/openfoam/-/blob/master/doc/Build.md

# Compiling uniGasFoam
**Source the openFOAM enviroment:**

Prior to compiling uniGasFoam source the correct OpenFOAM environment. For example, for the OpenFOAM-v2212 version:

* source [installation path]/OpenFOAM-v2212/etc/bashrc

**Compile libraries and executables**

From within the main repository directory execute the Allwamke script:
* ./Allwmake [nprocs]


# Running uniGasFoam
The tutorial cases found in the *tutorials/uniGasFoam* directory are a great starting point to run and test uniGasFoam.

To run uniGasFoam as a background process and output the case progress to a log file, execute:
* uniGasFoam > log &

Like any other OpenFOAM solver, uniGasFoam can be run in parallel by decomposing the flow domain. To decompose the case and run uniGasFoam in parallel in the background and output the case progress to a log file, execute in sequence:
* decomposePar
* uniGasFoam -n [nprocs] -parallel > log &
