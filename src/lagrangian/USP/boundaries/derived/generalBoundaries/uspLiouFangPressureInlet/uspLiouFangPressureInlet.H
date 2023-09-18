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
    Foam::uspLiouFangPressureInlet

Description
    This class provides functions to implement an implicit treatment for low
    speed inflow boundaries. The local mean flow velocity, temperature, and
    density are used to determine the number of molecules entering the
    computational domain and their corresponding velocities at every sample
    average step.

    Use in the boundariesDict file, in conjuction with a deletion patch.

    Reference: W.W. Liou and Y.C. Fang. "Implicit Boundary Conditions for Direct Simulation Monte Carlo Method in MEMS Flow Predictions", CMES 1(4), 110-128, 2000.

Usage
    \table
        Property         | Description           | Required  | Default
        typeIds          | List of species       | yes       |
        inletPressure    | Inlet pressure [Pa]   | yes       |
        inletTemperature | Inlet temperature [K] | yes       |
        moleFractions    | Mole fractions        | yes       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    generalBoundaryProperties
    {
        patchName   inlet;
    }

    boundaryModel   uspLiouFangPressureInlet;

    uspLiouFangPressureInletProperties
    {
        inletPressure       252000;
        inletTemperature    300;
        typeIds             (N2);
        moleFractions
        {
            N2      1.0;
        }
    }
    \endverbatim

SourceFiles
    uspLiouFangPressureInlet.C

\*---------------------------------------------------------------------------*/

#ifndef uspLiouFangPressureInlet_H
#define uspLiouFangPressureInlet_H

#include "uspGeneralBoundary.H"
#include "polyMesh.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class uspLiouFangPressureInlet Declaration
\*---------------------------------------------------------------------------*/

class uspLiouFangPressureInlet
:
    public uspGeneralBoundary
{
    // Private data

        dictionary propsDict_;
        scalarField moleFractions_;
        vectorField inletVelocity_;
        vectorField previousInletVelocity_;
        scalar inletPressure_;
        scalar inletTemperature_;
        scalar theta_;
        scalar n_;


    // Private Member Functions

        void setProperties();


public:

    //- Runtime type information
    TypeName("uspLiouFangPressureInlet");

    //- Construct from components
    uspLiouFangPressureInlet
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspLiouFangPressureInlet() = default;


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
