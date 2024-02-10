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
    Foam::uspFreeStreamInflowPatch

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

    boundaryModel   uspFreeStreamInflowPatch;

    uspFreeStreamInflowPatchProperties
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
    uspFreeStreamInflowPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uspFreeStreamInflowPatch_H
#define uspFreeStreamInflowPatch_H

#include "uspGeneralBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class uspFreeStreamInflowPatch Declaration
\*---------------------------------------------------------------------------*/

class uspFreeStreamInflowPatch
:
    public uspGeneralBoundary
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
    TypeName("uspFreeStreamInflowPatch");

    //- Construct from components
    uspFreeStreamInflowPatch
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspFreeStreamInflowPatch() = default;


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
