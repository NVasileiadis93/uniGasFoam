/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

Class
    uniGasZoneFill

Description
    This class provides functions to initialise particles based
    on uniform macroscopic properties.

    Use in the uniGasInitialisationDict file.

Usage
    Example of the particle initialization specification:
    
    \verbatim
    configuration
    {
        type                        uniGasZoneFill;
        zone                        uniGasZone;
        numberDensities{Ar          4.24700E+20;};
		translationalTemperature    200;
		rotationalTemperature       0;
		vibrationalTemperature      0;
        electronicTemperature       0;
		velocity		            (2634.7 0 0);
    }
    \endverbatim

SourceFiles
    uniGasZoneFill.C


\*---------------------------------------------------------------------------*/

#ifndef uniGasZoneFill_H
#define uniGasZoneFill_H

#include "uniGasConfiguration.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class uniGasZoneFill Declaration
\*---------------------------------------------------------------------------*/

class uniGasZoneFill
:
    public uniGasConfiguration
{
public:

    //- Runtime type information
    TypeName("uniGasZoneFill");

    //- Construct from components
    uniGasZoneFill(uniGasCloud& cloud, const dictionary& dict);

    //- Destructor
    virtual ~uniGasZoneFill() = default;


    // Member Functions

        virtual void setInitialConfiguration();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
