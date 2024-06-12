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
    Foam::noTimeCounterSubCycled

Description
    The no-time-counter collision partner selection routine, which is sub-cycled to account for the presence of electrons.

    Use in the uniGasProperties file.

    References:
        G.A. Bird,  "Molecular Gas Dynamics and the Direct Simulation of Gas Flows", Oxford Science Publications, 1994 (page 219).
        
        I.D. Boyd. "Computation of Hypersonic Flows Using the Direct Simulation Monte Carlo Method", Journal of Spacecraft and Rockets, 52(1), 38–53, 2015.
        
Usage
    Example of the binary collision partner selection specification:
    
    \verbatim
    binaryCollisionPartnerModel   		 noTimeCounterSubCycled;
    
    collisionProperties
    {
        nSubCycles      3;
    }
    \endverbatim

SourceFiles
    noTimeCounterSubCycled.C

\*---------------------------------------------------------------------------*/

#ifndef noTimeCounterSubCycled_H
#define noTimeCounterSubCycled_H


#include "binaryCollisionPartner.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class noTimeCounterSubCycled Declaration
\*---------------------------------------------------------------------------*/

class noTimeCounterSubCycled
:
    public binaryCollisionPartner
{
private:

    // Private data

        dictionary propsDict_;

        label nSubCycles_;

        label infoCounter_;


public:

    //- Runtime type information
    TypeName("noTimeCounterSubCycled");

    //- Construct from components
    noTimeCounterSubCycled
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uniGasCloud& cloud
    );

    //- Destructor
    virtual ~noTimeCounterSubCycled() = default;

    // Member Functions

        void initialConfiguration();

        void collide();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
