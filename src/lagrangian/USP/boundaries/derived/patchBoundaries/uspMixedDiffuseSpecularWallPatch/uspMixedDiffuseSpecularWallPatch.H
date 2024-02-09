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
    uspMixedDiffuseSpecularWallPatch

Description
    Wall interaction that performs a diffuse reflection or a specular
    reflection depending on the accommodation factor compared to a randomly
    generated number.

    Use in the boundariesDict file.

Usage

    Example of the boundary condition specification:
    \verbatim
    patchBoundaryProperties
    {
        patch               wall;
    }

    boundaryModel   uspMixedDiffuseSpecularWallPatch;

    uspMixedDiffuseSpecularWallPatchProperties
    {
        diffuseFraction     0.9;
        temperature         300;
        velocity            (0 0 0);
    }
    \endverbatim

SourceFiles
    uspMixedDiffuseSpecularWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uspMixedDiffuseSpecularWallPatch_H
#define uspMixedDiffuseSpecularWallPatch_H

#include "uspDiffuseWallPatch.H"
#include "uspSpecularWallPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class uspMixedDiffuseSpecularWallPatch Declaration
\*---------------------------------------------------------------------------*/

class uspMixedDiffuseSpecularWallPatch
:
    public uspPatchBoundary
{
    // Private data

        dictionary propsDict_;

        //- Fraction of wall interactions that are diffuse
        scalar diffuseFraction_;

        scalar temperature_;

        vector velocity_;

public:

    //- Runtime type information
    TypeName("uspMixedDiffuseSpecularWallPatch");

    // Constructors

        //- Construct from components
        uspMixedDiffuseSpecularWallPatch
        (
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        );

    //- Destructor
    virtual  ~uspMixedDiffuseSpecularWallPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control the mols
        virtual void controlParticle(uspParcel& p, uspParcel::trackingData& td);

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
