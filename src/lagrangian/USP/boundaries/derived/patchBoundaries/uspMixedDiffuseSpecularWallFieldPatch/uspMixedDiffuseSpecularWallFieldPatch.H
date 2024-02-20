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
    Foam::uspMixedDiffuseSpecularWallFieldPatch

Description
    Wall interaction that performs a diffuse reflection or a specular
    reflection depending on the accommodation factor compared to a randomly
    generated number.

    Use in the boundariesDict file.

    FIELD VERSION: require boundaryU and boundaryT fields to be defined in
    your time directory.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    patchBoundaryProperties
    {
        patch               wall;
    }

    boundaryModel   uspMixedDiffuseSpecularWallFieldPatch;

    uspMixedDiffuseSpecularWallFieldPatchProperties
    {
        diffuseFraction     0.9;
    }
    \endverbatim

SourceFiles
    uspMixedDiffuseSpecularWallFieldPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uspMixedDiffuseSpecularWallFieldPatch_H
#define uspMixedDiffuseSpecularWallFieldPatch_H

#include "uspDiffuseWallPatch.H"
#include "uspSpecularWallPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
           Class uspMixedDiffuseSpecularWallFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class uspMixedDiffuseSpecularWallFieldPatch
:
    public virtual uspPatchBoundary
{
    // Private data

        dictionary propsDict_;

        //- Fraction of wall interactions that are diffuse
        scalar diffuseFraction_;

        //- Boundary temperature
        volScalarField boundaryT_;

        //- Boundary velocity
        volVectorField boundaryU_;

public:

    //- Runtime type information
    TypeName("uspMixedDiffuseSpecularWallFieldPatch");

    // Constructors

        //- Construct from components
        uspMixedDiffuseSpecularWallFieldPatch
        (
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        );

    //- Destructor
    virtual ~uspMixedDiffuseSpecularWallFieldPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- calculate any required properties
        virtual void calculateProperties();

        //- control the mols
        virtual void controlParticle(uspParcel& p, uspParcel::trackingData& td);

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
