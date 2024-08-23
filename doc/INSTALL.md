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
