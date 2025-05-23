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
    Foam::noTimeCounter

Description
    The no-time-counter collision partner selection routine.

    Use in the uniGasProperties file.

    Reference:
        G.A. Bird,  "Molecular Gas Dynamics and the Direct Simulation of Gas Flows", Oxford Science Publications, 1994 (page 219).
        
Usage
    Example of the dsmc collision partner selection specification:

    \verbatim
    dsmcCollisionPartnerModel   		 noTimeCounter;
    \endverbatim

SourceFiles
    noTimeCounter.C

\*---------------------------------------------------------------------------*/

#ifndef noTimeCounter_H
#define noTimeCounter_H

#include "dsmcCollisionPartner.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class noTimeCounter Declaration
\*---------------------------------------------------------------------------*/

class noTimeCounter
:
    public dsmcCollisionPartner
{
    // Private data

        label infoCounter_;


public:

    //- Runtime type information
    TypeName("noTimeCounter");

    //- Construct from components
    noTimeCounter
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uniGasCloud& cloud
    );

    //- Destructor
    virtual ~noTimeCounter() = default;

    // Member Functions

        void initialConfiguration();

        void collide();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
