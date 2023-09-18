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
    uspCLLWallFieldPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    uspCLLWallFieldPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uspCLLWallFieldPatch_H
#define uspCLLWallFieldPatch_H

#include "uspPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "uspCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class uspCLLWallFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class uspCLLWallFieldPatch
:
    public uspPatchBoundary
{
    // Private data

        dictionary propsDict_;

        void setProperties();

        scalar normalAccommodationCoefficient_;

        scalar tangentialAccommodationCoefficient_;

        scalar rotationalEnergyAccommodationCoefficient_;

        //- Boundary temperature
        volScalarField boundaryT_;

        //- Boundary velocity
        volVectorField boundaryU_;


public:

    //- Runtime type information
    TypeName("uspCLLWallFieldPatch");

    //- Construct from components
    uspCLLWallFieldPatch
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspCLLWallFieldPatch() = default;


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
