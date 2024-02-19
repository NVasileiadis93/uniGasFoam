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
    Foam::uspMassFlowRateInletPatch

Description
    This class provides functions to implement an a known mass flow rate for 
    low speed outflow boundaries. The local mean flow velocity, temperature, and
    density are used to determine the number of molecules entering the
    computational domain and their corresponding velocities at every sample
    average step.

    Use in the boundariesDict file, in conjunction with a deletion patch.

    Reference:
        M. Lei et al,  "The implementation of subsonic boundary conditions 
        for the direct simulation Monte Carlo method in uspFoam", 
        Computers & Fluids, 156, 209–219, 2017.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    generalBoundaryProperties
    {
        patch       inlet;
    }

    boundaryModel   uspMassFlowRateInletPatch;

    uspMassFlowRateInletPatchProperties
    {
        typeIds                     (Ar);
        moleFractions               {Ar     1.0;}
        inletTemperature            300;
        massFlowRate                1e-11;
        initialVelocity             (1 0 0);
        theta                       0.01;
    }
    \endverbatim

SourceFiles
    uspMassFlowRateInletPatch.C

\*---------------------------------------------------------------------------*/

#ifndef uspMassFlowRateInletPatch_H
#define uspMassFlowRateInletPatch_H

#include "uspGeneralBoundary.H"
#include "polyMesh.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class uspMassFlowRateInletPatch Declaration
\*---------------------------------------------------------------------------*/

class uspMassFlowRateInletPatch
:
    public uspGeneralBoundary
{
    // Private data

        dictionary propsDict_;

        scalarField moleFractions_;
        
        scalar massFlowRate_;
        
        scalar inletTemperature_;
        
        vector initialVelocity_;
        
        scalar theta_;

        scalarField moleFlowRate_;
        
        scalarField parcelsIn_;
        
        scalarField parcelsOut_;
        
        scalarField parcelsToInsert_;
        
        List<scalarField> inletNumberDensity_;
        
        vectorField inletVelocity_;
        
        vectorField previousInletVelocity_;

public:

    //- Runtime type information
    TypeName("uspMassFlowRateInletPatch");

    //- Construct from components
    uspMassFlowRateInletPatch
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspMassFlowRateInletPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control the mols
        virtual void controlParcelsBeforeMove();
        virtual void controlParcelsBeforeCollisions();
        virtual void controlParcelsAfterCollisions();

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
