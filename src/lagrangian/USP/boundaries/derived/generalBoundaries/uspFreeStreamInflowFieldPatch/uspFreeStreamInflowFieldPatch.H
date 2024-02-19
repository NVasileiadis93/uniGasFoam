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
    Foam::uspFreeStreamInflowFieldPatch

Description
    This class provides functions to insert the correct number density,
    velocity, and temperature at a freestream inlet characterized by the
    Maxwell distribtution. 
    
    Use in the boundariesDict file, in conjuction with a deletion patch.

    FIELD VERSION: require numberDensity_Ar, boundaryTransT, boundaryRotT, 
    boundaryVibT, boundaryElecT and boundaryU fields to be defined in
    your time directory.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    generalBoundaryProperties
    {
        patch       flow;
    }

    boundaryModel   uspFreeStreamInflowFieldPatch;

    uspFreeStreamInflowFieldPatchProperties
    {
        typeIds     (Ar);
    }
    \endverbatim

SourceFiles
    uspFreeStreamInflowFieldPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uspFreeStreamInflowFieldPatch_H
#define uspFreeStreamInflowFieldPatch_H

#include "uspFreeStreamInflowPatch.H"
#include "polyMesh.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class uspFreeStreamInflowFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class uspFreeStreamInflowFieldPatch
:
    public uspGeneralBoundary
{
    // Private data

        dictionary propsDict_;

        // List of (one value) mass densities corresponding to the type Ids
        List<scalarField> numberDensities_;

        // Translational temperature
        scalarField inletTransT_;

        // Rotational temperature
        scalarField inletRotT_;

        // Vibrational temperature
        scalarField inletVibT_;

        // Electronic temperature
        scalarField inletElecT_;

        vectorField inletVelocities_;

        //- List of boundary number densities (one entry for each species)
        List<autoPtr<volScalarField>> boundaryNumberDensity_;

        //- Boundary temperatures
        volScalarField boundaryTransT_;

        volScalarField boundaryRotT_;

        volScalarField boundaryVibT_;

        volScalarField boundaryElecT_;

        //- Boundary velocity
        volVectorField boundaryU_;

public:

    //- Runtime type information
    TypeName("uspFreeStreamInflowFieldPatch");

    //- Construct from components
    uspFreeStreamInflowFieldPatch
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    // Destructor
    virtual ~uspFreeStreamInflowFieldPatch() = default;


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
