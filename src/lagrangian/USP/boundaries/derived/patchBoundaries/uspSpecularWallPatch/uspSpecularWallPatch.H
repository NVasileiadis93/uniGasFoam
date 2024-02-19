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
    Foam::uspSpecularWallPatch

Description

     Specular wall boundary patch: reflecting molecules back in by reversing
     their normal velocity vector to the boundary.

     Use in the boundariesDict file.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    patchBoundaryProperties
    {
        patchName                           wall;
    }

    boundaryModel   uspSpecularWallPatch;

    uspSpecularWallPatchProperties
    {
    }
    \endverbatim

SourceFiles
    uspSpecularWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uspSpecularWallPatch_H
#define uspSpecularWallPatch_H

#include "uspPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class uspSpecularWallPatch Declaration
\*---------------------------------------------------------------------------*/

class uspSpecularWallPatch
:
    public virtual uspPatchBoundary
{
    // Private data

        dictionary propsDict_;


public:

    //- Runtime type information
    TypeName("uspSpecularWallPatch");

    // Constructors

        //- Construct from components
        uspSpecularWallPatch
        (
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        );

    //- Destructor
    virtual ~uspSpecularWallPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control the mols
        virtual void controlParticle
        (
            uspParcel& p,
            uspParcel::trackingData& td
        );

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
