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
    Foam::uniGasChapmanEnskogFreeStreamInflowPatch

Description
    This class provides functions to insert the correct number density,
    velocity, and temperature at a freestream inlet characterized by the
    Chapman-Enskog distribtution. This class is for uniform properties 
    across the patch.

    Use in the boundariesDict file, in conjuction with a deletion patch.

Usage
    Example of the boundary condition specification:
    
    \verbatim
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
    \endverbatim

SourceFiles
    uniGasChapmanEnskogFreeStreamInflowPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uniGasChapmanEnskogFreeStreamInflowPatch_H
#define uniGasChapmanEnskogFreeStreamInflowPatch_H

#include "uniGasGeneralBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class uniGasChapmanEnskogFreeStreamInflowPatch Declaration
\*---------------------------------------------------------------------------*/

class uniGasChapmanEnskogFreeStreamInflowPatch
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

        vector heatFlux_;
    
        tensor stress_;

public:

    //- Runtime type information
    TypeName("uniGasChapmanEnskogFreeStreamInflowPatch");

    //- Construct from components
    uniGasChapmanEnskogFreeStreamInflowPatch
    (
        const polyMesh& mesh,
        uniGasCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uniGasChapmanEnskogFreeStreamInflowPatch() = default;


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