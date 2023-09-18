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
    Foam::uspFreeStreamInflowPatch

Description
    This class provides functions to insert the correct number density,
    velocity, and temperature at a freestream inlet. This class is
    for uniform properties across the patch.

    From Bird's 1994 book, equation 4.22

        \f[
            \dot{N}_{i} = \left(n c_{mp}\right)\frac{[\left[ exp(-s^2\cos^2\theta + \pi^(1/2)\cos\theta\left{ 1 + \ref\left(s\cos\theta \right)\right}\right]}{2\pi^(1/2)}
        \f]

    where
    \vartable
        n           | Number density [#/m^3]
        c_{mp}      | Most probable thermal velocity [m/s]
        s           | Molecular speed ratio
        \theta      | Angle between face normal and inlet velocity vectors
        \dot{N}_{i} | Inlet number flux [#/m^2/s]
    \endvartable

    Use in the boundariesDict file, in conjuction with a deletion patch.

    Reference:
        G.A. Bird,  "Molecular Gas Dynamics and the Direct Simulation of Gas Flows", Oxford Science Publications, 1994 (page 82).

Usage
    \table
        Property                 | Description                   | Required  | Default
        typeIds                  | List of species               | yes       |
        translationalTemperature | Translational temperature [K] | yes       |
        rotationalTemperature    | Rotational temperature [K]    | yes       |
        vibrationalTemperature   | Vibrational temperature [K]   | yes       |
        electronicTemperature    | Electronic temperature [K]    | yes       |
        velocity                 | Inlet velocity [m/s]          | yes       |
        numberDensities          | Number density list [#/m^3]   | yes       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    generalBoundaryProperties
    {
        patchName                           flow;
    }

    boundaryModel   uspFreeStreamInflowPatch;

    uspFreeStreamInflowPatchProperties
    {
        typeIds                     (Ar);
        translationalTemperature    100;
        rotationalTemperature       100;
        vibrationalTemperature      100;
        electronicTemperature       100;
        velocity                    (1000 0 0);
        numberDensities
        {
            Ar          1.0e21;
        }
    }
    \endverbatim

SourceFiles
    uspFreeStreamInflowPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uspFreeStreamInflowPatch_H
#define uspFreeStreamInflowPatch_H

#include "uspGeneralBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class uspFreeStreamInflowPatch Declaration
\*---------------------------------------------------------------------------*/

class uspFreeStreamInflowPatch
:
    public uspGeneralBoundary
{
    // Private data

        dictionary propsDict_;
        void setProperties();


public:

    //- Runtime type information
    TypeName("uspFreeStreamInflowPatch");

    //- Construct from components
    uspFreeStreamInflowPatch
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspFreeStreamInflowPatch() = default;


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
