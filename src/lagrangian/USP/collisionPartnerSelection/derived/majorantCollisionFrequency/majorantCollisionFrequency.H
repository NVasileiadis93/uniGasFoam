/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    majorantCollisionFrequency

Description

The majorant collision frequency method. 
See https://aip.scitation.org/doi/10.1063/1.5119580 for details.

SourceFiles
    majorantCollisionFrequencyI.H
    majorantCollisionFrequency.C
    majorantCollisionFrequencyIO.C

\*---------------------------------------------------------------------------*/

#ifndef majorantCollisionFrequency_H
#define majorantCollisionFrequency_H


#include "collisionPartnerSelection.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class majorantCollisionFrequency Declaration
\*---------------------------------------------------------------------------*/

class majorantCollisionFrequency
:
    public collisionPartnerSelection
{

private:

    // Private data

        label infoCounter_;

public:

    //- Runtime type information
    TypeName("majorantCollisionFrequency");

    // Static data members

    // Constructors

        //- construct from components
        majorantCollisionFrequency
        (
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~majorantCollisionFrequency();

    // Member Functions

        // Access

        // Check

        // Edit

            void initialConfiguration();

            void collide();

        // Write





    // Member Operators

//         void operator=(const majorantCollisionFrequency&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const majorantCollisionFrequency&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
