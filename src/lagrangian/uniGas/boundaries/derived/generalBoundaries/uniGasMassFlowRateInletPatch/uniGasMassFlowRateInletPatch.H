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
    Foam::uniGasMassFlowRateInletPatch

Description
    This class provides functions to implement an a known mass flow rate for 
    low speed outflow boundaries. The local mean flow velocity, temperature, and
    density are used to determine the number of molecules entering the
    computational domain and their corresponding velocities at every sample
    average step.

    Use in the boundariesDict file, in conjunction with a deletion patch.

    Reference:
        M. Lei et al,  "The implementation of subsonic boundary conditions 
        for the direct simulation Monte Carlo method in uniGasFoam", 
        Computers & Fluids, 156, 209–219, 2017.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    generalBoundaryProperties
    {
        patch       inlet;
    }

    boundaryModel   uniGasMassFlowRateInletPatch;

    uniGasMassFlowRateInletPatchProperties
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
    uniGasMassFlowRateInletPatch.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasMassFlowRateInletPatch_H
#define uniGasMassFlowRateInletPatch_H

#include "uniGasGeneralBoundary.H"
#include "polyMesh.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class uniGasMassFlowRateInletPatch Declaration
\*---------------------------------------------------------------------------*/

class uniGasMassFlowRateInletPatch
:
    public uniGasGeneralBoundary
{
    // Private data

        dictionary propsDict_;
        scalar theta_;
        scalar inletTemperature_;
        scalar massFlowRate_;
        vector initialVelocity_;
        scalarField moleFractions_;
        scalar patchSurfaceArea_;
        scalarField moleFlowRate_;
        scalarField parcelsIn_;
        scalarField parcelsToInsert_;
        List<scalarField> inletNumberDensity_;
        vectorField inletVelocity_;
        vectorField previousInletVelocity_;


    // Private Member Functions

        void setProperties();


public:

    //- Runtime type information
    TypeName("uniGasMassFlowRateInletPatch");

    //- Construct from components
    uniGasMassFlowRateInletPatch
    (
        const polyMesh& mesh,
        uniGasCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uniGasMassFlowRateInletPatch() = default;


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