/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Foam::noTimeCounterSubCycled

Description
    The no-time-counter collision partner selection routine, which is sub-cycled to account for the presence of electrons.
    
    From Bird's 1994 book, equation 11.3, the NTC method is that

        \f[
            0.5 N (N-1) F_N \left( (\sigma_T c_r\right)_{max} \Delta t/V_C
        \f]

    pairs are selected from a cell each time step, where
    \vartable
        N           | Number of USP particles in the cell
        F_N         | Number of real atoms/molecules the USP particles represent
        \sigma_T    | Collision cross-section [m^2]
        c_r         | Relative speed of the collision pair [m/s]
        \Delta t    | Time step [s]
        V_C         | Cell volume [m^3]
    \endvartable
    
    The pair are selected for collision if
    
        \f[
            \frac{\sigma_T c_r}{\left(\sigma_T c_r\right)_{max} > Rand}
        \f]
        
    where
    \vartable
        Rand         | Uniformly distributed randomly generated number between 0 and 1.
    \endvartable
    
    References:
        G.A. Bird,  "Molecular Gas Dynamics and the Direct Simulation of Gas Flows", Oxford Science Publications, 1994 (page 219).
        
        I.D. Boyd. "Computation of Hypersonic Flows Using the Direct Simulation Monte Carlo Method", Journal of Spacecraft and Rockets, 52(1), 38–53, 2015.
        
Usage
    \verbatim
    binaryCollisionPartnerModel   		 noTimeCounterSubCycled;
    
    noTimeCounterSubCycledProperties
    {
        nSubCycles      30;
    }
    \endverbatim

SourceFiles
    noTimeCounterSubCycledI.H
    noTimeCounterSubCycled.C

\*---------------------------------------------------------------------------*/

#ifndef noTimeCounterSubCycled_H
#define noTimeCounterSubCycled_H


#include "binaryCollisionPartner.H"
#include "uspCloud.H"

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
        uspCloud& cloud
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
