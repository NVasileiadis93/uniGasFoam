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
    Foam::uspDiffuseWallFieldPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature and velocity

    Use in the boundariesDict file.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    patchBoundaryProperties
    {
        patch           wall;
    }

    boundaryModel   uspDiffuseWallFieldPatch;

    uspDiffuseWallFieldPatchProperties
    {
    }
    \endverbatim

    FIELD VERSION: require boundaryU and boundaryT fields to be defined in
    time directory.

SourceFiles
    uspDiffuseWallFieldPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uspDiffuseWallFieldPatch_H
#define uspDiffuseWallFieldPatch_H

#include "uspDiffuseWallPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class uspDiffuseWallFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class uspDiffuseWallFieldPatch
:
    public uspPatchBoundary
{
    // Private data

        dictionary propsDict_;

        //- Boundary temperature
        volScalarField boundaryT_;

        //- Boundary velocity
        volVectorField boundaryU_;


public:

    //- Runtime type information
    TypeName("uspDiffuseWallFieldPatch");

    // Constructors

        //- Construct from components
        uspDiffuseWallFieldPatch
        (
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        );

    //- Destructor
    virtual ~uspDiffuseWallFieldPatch() = default;


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
