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
    Foam::uspLiouFangPressureOutletPatch

Description
    This class provides functions to implement an implicit treatment for low
    speed outflow boundaries. The local mean flow velocity, temperature, and
    density are used to determine the number of molecules entering the
    computational domain and their corresponding velocities at every sample
    average step.

    Use in the boundariesDict file, in conjuction with a deletion patch.

    Reference: W.W. Liou and Y.C. Fang. "Implicit Boundary Conditions for Direct
    Simulation Monte Carlo Method in MEMS Flow Predictions", CMES 1(4),
    110-128, 2000.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    generalBoundaryProperties
    {
        patch       outlet;
    }

    boundaryModel   uspLiouFangPressureOutletPatch;

    uspLiouFangPressureOutletPatchProperties
    {
        typeIds             (Ar);
        moleFractions       {Ar     1.0;}
        outletPressure      100000;
    }
    \endverbatim

SourceFiles
    uspLiouFangPressureOutletPatch.C

\*---------------------------------------------------------------------------*/

#ifndef uspLiouFangPressureOutletPatch_H
#define uspLiouFangPressureOutletPatch_H

#include "uspGeneralBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                 Class uspLiouFangPressureOutletPatch Declaration
\*---------------------------------------------------------------------------*/

class uspLiouFangPressureOutletPatch
:
    public uspGeneralBoundary
{
    // Private data

        dictionary propsDict_;
        
        scalar outletPressure_;

        label nTimeSteps_;
        vectorField UMean_;
        vectorField outletVelocity_;
        vectorField UCollected_;
        List<scalarField> totalVibrationalEnergy_;
        List<scalarField> vibT_;
        List<scalarField> vDof_;
        List<scalarField> nTotalParcelsSpecies_;
        scalarField moleFractions_;
        scalarField cellVolume_;
        scalarField outletNumberDensity_;
        scalarField outletMassDensity_;
        scalarField outletTemperature_;
        vectorField totalMomentum_;
        scalarField totalMass_;
        scalarField totalMassDensity_;
        scalarField totalNumberDensity_;
        scalarField totalRotationalEnergy_;
        scalarField totalRotationalDof_;
        scalarField nTotalParcels_;
        scalarField nTotalParcelsInt_;
        vectorField velSqrMean_;
        vectorField velMeanSqr_;
        vectorField totalVelSqrMean_;
        vectorField totalVelMeanSqr_;

public:

    //- Runtime type information
    TypeName("uspLiouFangPressureOutletPatch");

    //- Construct from components
    uspLiouFangPressureOutletPatch
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspLiouFangPressureOutletPatch() = default;


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
