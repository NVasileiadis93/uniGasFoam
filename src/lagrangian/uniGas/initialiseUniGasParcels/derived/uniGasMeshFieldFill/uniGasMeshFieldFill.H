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

Class
    Foam::uniGasMeshFieldFill

Description
    This class provides functions to initialise particles.

    Use in the uniGasInitialiseDict file.

    FIELD VERSION: require initialNumberDensity_Ar, initialTransT, initialRotT, 
    initialVibT, initialElecT and initialU fields to be defined in your time directory.

Usage
    Example of the particle initialization specification:
    
    \verbatim
    configuration
    {
        type            uniGasMeshFieldFill;
        typeIdList      (Ar);
    }
    \endverbatim

SourceFiles
    uniGasMeshFieldFill.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasMeshFieldFill_H
#define uniGasMeshFieldFill_H

#include "uniGasAllConfigurations.H"
#include "vector.H"
#include "polyMesh.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class uniGasMeshFieldFill Declaration
\*---------------------------------------------------------------------------*/

class uniGasMeshFieldFill
:
    public uniGasConfiguration
{
public:

    //- Runtime type information
    TypeName("uniGasMeshFieldFill");

    //- Construct from components
    uniGasMeshFieldFill(uniGasCloud& cloud, const dictionary& dict);

    //- Destructor
    virtual ~uniGasMeshFieldFill() = default;


    // Member Functions
        virtual void setInitialConfiguration();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
