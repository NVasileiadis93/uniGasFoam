# About uniGasFoam
The uniGasFoam solver is an open-source particle-based simulation framework for multiscale rarefied gas flows. It is entirely designed in OpenFOAM and is an extension of the well-established direct simulation Monte Carlo (DSMC) solver dsmcFoam+ (https://github.com/MicroNanoFlows/OpenFOAM-2.4.0-MNF/tree/master).

# Copyright
uniGasFoam is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See *doc/LICENSE* or http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you can copy the files.

# Installation and usage
Detailed instructions for installing and building uniGasFoam can be found in *doc/INSTALL.md*. Using uniGasFoam will be relatively straightforward for users who are already familiar with the OpenFOAM suite. A detailed description of case file structure and dictionary entries is provided in *doc/GUIDE.md*. as a guide to the uniGasFoam solver. 

The source files for the libraries are located in *src/lagrangian* and the executable file for running the solver can be found in *application/solvers/discreteMethods/uniGasFoam*. In addition, tutorial cases can be found in the *tutorials/uniGasFoam* directory.

# Further information
For bug reporting, code development and any further questions, feel free to reach me at [Nikos Vasieliadis](mailto:nikovasi93@gmail.com).
