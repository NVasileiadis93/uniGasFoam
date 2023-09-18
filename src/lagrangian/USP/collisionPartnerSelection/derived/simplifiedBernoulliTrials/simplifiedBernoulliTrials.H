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
    simplifiedBernoulliTrials

Description

Stefan Stefanov's simplified Bernoulli trials collision partner selection routine.
Leads to fewer repeat collisions than Bird's no time counter method.
See SIAM J. Sci. Comput. 33(2) 667-702 for details.
We do not use the half timestep implementation.

SourceFiles
    simplifiedBernoulliTrials.C

\*---------------------------------------------------------------------------*/

#ifndef simplifiedBernoulliTrials_H
#define simplifiedBernoulliTrials_H


#include "collisionPartnerSelection.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simplifiedBernoulliTrials Declaration
\*---------------------------------------------------------------------------*/

class simplifiedBernoulliTrials
:
    public collisionPartnerSelection
{

private:

    // Private data

        label infoCounter_;

        List<List<scalar>> subcellVolume;

public:

    //- Runtime type information
    TypeName("simplifiedBernoulliTrials");

    // Static data members

    // Constructors

        //- construct from components
        simplifiedBernoulliTrials
        (
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        );

    // Selectors


    // Destructor

        ~simplifiedBernoulliTrials();

    // Member Functions

        // Access

        // Check

        // Edit

            void initialConfiguration();

            void collide();

        // Write





    // Member Operators

//         void operator=(const simplifiedBernoulliTrials&);

    // Friend Functions

    // Friend Operators

    // IOstream Operators

//         friend Ostream& operator<<(Ostream&, const simplifiedBernoulliTrials&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
