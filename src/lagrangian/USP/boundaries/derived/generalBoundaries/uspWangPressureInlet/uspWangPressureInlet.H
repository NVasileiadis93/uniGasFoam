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
    Foam::uspWangPressureInlet

Description
    This class provides functions to implement an implicit treatment for low
    speed outflow boundaries. The local mean flow velocity, temperature, and
    density are used to determine the number of molecules entering the
    computational domain and their corresponding velocities at every sample
    average step.

    Use in the boundariesDict file, in conjuction with a deletion patch.

    Reference: M. Wang and Z. Li. "Simulations for gas flows in microgeometries using the direct simulation Monte Carlo method". International Journal of Heat and Fluid Flow, 25(6):975–985, 2004

Usage
    \table
        Property         | Description           | Required  | Default
        typeIds          | List of species       | yes       |
        outletPressure   | outlet pressure [Pa]  | yes       |
        moleFractions    | Mole fractions        | yes       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    generalBoundaryProperties
    {
        patchName   inlet;
    }

    boundaryModel   uspWangPressureInlet;

    uspWangPressureInletProperties
    {
        outletPressure       252000;
        typeIds             (N2);
        moleFractions
        {
            N2      1.0;
        }
    }
    \endverbatim

SourceFiles
    uspWangPressureInlet.C


\*---------------------------------------------------------------------------*/

#ifndef uspWangPressureInlet_H
#define uspWangPressureInlet_H

#include "uspGeneralBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class uspWangPressureInlet Declaration
\*---------------------------------------------------------------------------*/

class uspWangPressureInlet
:
    public uspGeneralBoundary
{
    // Private data

        dictionary propsDict_;

        scalarField moleFractions_;
        scalar inletPressure_;
        scalar inletTemperature_;
        scalar n_;
        scalarField cellVolume_;
        vectorField inletVelocity_;
        vectorField totalMomentum_;
        scalarField totalMass_;
        scalarField nTotalParcels_;
        scalar nTimeSteps_;
        scalarField mcc_;
        vectorField UMean_;
        vectorField UCollected_;
        vectorField velSqrMean_;
        vectorField velMeanSqr_;
        vectorField totalVelSqrMean_;
        vectorField totalVelMeanSqr_;


    // Private Member Functions

        void setProperties();


public:

    //- Runtime type information
    TypeName("uspWangPressureInlet");

    //- Construct from components
    uspWangPressureInlet
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspWangPressureInlet() = default;


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
