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

Class
    Foam::uspMeshFieldFill

Description
    This class provides functions to initialise particles.

    Use in the uspInitialiseDict file.

    FIELD VERSION: require initialNumberDensity_Ar, initialTransT, initialRotT, 
    initialVibT, initialElecT and initialU fields to be defined in your time directory.

Usage
    Example of the particle initialization specification:
    
    \verbatim
    configuration
    {
        type            uspMeshFieldFill;
        typeIdList      (Ar);
    }
    \endverbatim

SourceFiles
    uspMeshFieldFill.C

\*---------------------------------------------------------------------------*/

#ifndef uspMeshFieldFill_H
#define uspMeshFieldFill_H

#include "uspAllConfigurations.H"
#include "vector.H"
#include "polyMesh.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class uspMeshFieldFill Declaration
\*---------------------------------------------------------------------------*/

class uspMeshFieldFill
:
    public uspConfiguration
{
public:

    //- Runtime type information
    TypeName("uspMeshFieldFill");

    //- Construct from components
    uspMeshFieldFill(uspCloud& cloud, const dictionary& dict);

    //- Destructor
    virtual ~uspMeshFieldFill() = default;


    // Member Functions
        virtual void setInitialConfiguration();
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
