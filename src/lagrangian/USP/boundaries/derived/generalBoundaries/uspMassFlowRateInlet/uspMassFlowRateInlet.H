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
    Foam::uspMassFlowRateInlet

Description
    This class provides functions to insert the mass flow rate at a patch,
    given a temperature

    Use in the boundariesDict file, in conjunction with a deletion patch.

    Reference:
        M. Lei et al,  "The implementation of subsonic boundary conditions for the direct simulation Monte Carlo method in uspFoam", Computers & Fluids, 156, 209–219, 2017.

    \table
        Property         | Description                        | Required  | Default
        theta            | Instanteneous velocity weight      | no        | 1.0
        typeIds          | List of species                    | yes       |
        inletTemperature | Inlet temperature [K]              | yes       |
        massFlowRate     | Mass flow rate [kg/s]              | yes       |
        initialVelocity  | Guess for the inlet velocity [m/s] | yes       |
        moleFractions    | Mole fractions                     | yes       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    generalBoundaryProperties
    {
        patchName                           flow;
    }

    boundaryModel   uspMassFlowRateInletPatch;

    uspMassFlowRateInletPatchProperties
    {
        theta                       0.01;
        typeIds                     (Ar);
        inletTemperature            300;
        massFlowRate                1e-11;
        initialVelocity             (1 0 0);
        moleFractions
        {
            Ar          1.0;
        }
    }
    \endverbatim

SourceFiles
    uspMassFlowRateInlet.C

\*---------------------------------------------------------------------------*/

#ifndef uspMassFlowRateInlet_H
#define uspMassFlowRateInlet_H

#include "uspGeneralBoundary.H"
#include "polyMesh.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class uspMassFlowRateInlet Declaration
\*---------------------------------------------------------------------------*/

class uspMassFlowRateInlet
:
    public uspGeneralBoundary
{
    // Private data

        dictionary propsDict_;
        scalar theta_;
        scalar inletTemperature_;
        scalar massFlowRate_;
        scalarField moleFractions_;
        scalarField moleFlowRate_;
        scalarField parcelsIn_;
        scalarField parcelsOut_;
        scalarField parcelsToInsert_;
        List<scalarField> n_;
        vector initialVelocity_;
        vectorField inletVelocity_;
        vectorField previousInletVelocity_;


    // Private Member Functions

        void setProperties();


public:

    //- Runtime type information
    TypeName("uspMassFlowRateInlet");

    //- Construct from components
    uspMassFlowRateInlet
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspMassFlowRateInlet() = default;


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
