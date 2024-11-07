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
    Foam::uniGasMixedDiffuseSpecularWallPatch

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

    boundaryModel   uniGasMixedDiffuseSpecularWallPatch;

    uniGasMixedDiffuseSpecularWallPatchProperties
    {
        diffuseFraction     0.9;
        temperature         300;
        velocity            (0 0 0);
    }
    \endverbatim

SourceFiles
    uniGasMixedDiffuseSpecularWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uniGasMixedDiffuseSpecularWallPatch_H
#define uniGasMixedDiffuseSpecularWallPatch_H

#include "uniGasDiffuseWallPatch.H"
#include "uniGasSpecularWallPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class uniGasMixedDiffuseSpecularWallPatch Declaration
\*---------------------------------------------------------------------------*/

class uniGasMixedDiffuseSpecularWallPatch
:
    public uniGasPatchBoundary
{
    // Private data

        dictionary propsDict_;

        //- Fraction of wall interactions that are diffuse
        scalar diffuseFraction_;

        scalar temperature_;

        vector velocity_;

public:

    //- Runtime type information
    TypeName("uniGasMixedDiffuseSpecularWallPatch");

    // Constructors

        //- Construct from components
        uniGasMixedDiffuseSpecularWallPatch
        (
            const polyMesh& mesh,
            uniGasCloud& cloud,
            const dictionary& dict
        );

    //- Destructor
    virtual  ~uniGasMixedDiffuseSpecularWallPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control the mols
        virtual void controlParticle(uniGasParcel& p, uniGasParcel::trackingData& td);

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
