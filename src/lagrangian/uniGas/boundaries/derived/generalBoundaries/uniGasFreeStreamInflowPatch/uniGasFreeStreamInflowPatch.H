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
    Foam::uniGasFreeStreamInflowPatch

Description
    This class provides functions to insert the correct number density,
    velocity, and temperature at a freestream inlet characterized by the
    Maxwell distribtution. This class is for uniform properties 
    across the patch.

    Use in the boundariesDict file, in conjuction with a deletion patch.

Usage
    Example of the boundary condition specification:
    
    \verbatim
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
    \endverbatim

SourceFiles
    uniGasFreeStreamInflowPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uniGasFreeStreamInflowPatch_H
#define uniGasFreeStreamInflowPatch_H

#include "uniGasGeneralBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class uniGasFreeStreamInflowPatch Declaration
\*---------------------------------------------------------------------------*/

class uniGasFreeStreamInflowPatch
:
    public uniGasGeneralBoundary
{
    // Private data

        dictionary propsDict_;

        scalar translationalTemperature_;
        
        scalar rotationalTemperature_;
    
        scalar vibrationalTemperature_;
    
        scalar electronicTemperature_;
    
        vector velocity_;

public:

    //- Runtime type information
    TypeName("uniGasFreeStreamInflowPatch");

    //- Construct from components
    uniGasFreeStreamInflowPatch
    (
        const polyMesh& mesh,
        uniGasCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uniGasFreeStreamInflowPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control the mols
        virtual void controlParcelsBeforeMove();
        virtual void controlParcelsBeforeCollisions();
        virtual void controlParcelsAfterCollisions();

        //- Output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- Update properties from a modified dictionary
        virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
