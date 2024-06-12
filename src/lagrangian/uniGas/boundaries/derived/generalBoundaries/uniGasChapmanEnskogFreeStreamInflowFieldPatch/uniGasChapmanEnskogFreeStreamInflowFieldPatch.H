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
    Foam::uniGasChapmanEnskogFreeStreamInflowFieldPatch

Description
    This class provides functions to insert the correct number density,
    velocity, and temperature at a freestream inlet characterized by the
    Chapman-Enskog distribtution.

    Use in the boundariesDict file, in conjuction with a deletion patch.

    FIELD VERSION: require numberDensity_Ar, boundaryTransT, boundaryRotT, 
    boundaryVibT, boundaryElecT, boundaryU, boundaryHeatFlux and
    boundaryStress fields to be defined in your time directory.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    generalBoundaryProperties
    {
        patch       flow;
    }

    boundaryModel   uniGasChapmanEnskogFreeStreamInflowFieldPatch;

    uniGasChapmanEnskogFreeStreamInflowFieldPatchProperties
    {
        typeIds     (Ar);
    }
    \endverbatim

SourceFiles
    uniGasChapmanEnskogFreeStreamInflowFieldPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uniGasChapmanEnskogFreeStreamInflowFieldPatch_H
#define uniGasChapmanEnskogFreeStreamInflowFieldPatch_H

#include "uniGasGeneralBoundary.H"
#include "uniGasChapmanEnskogFreeStreamInflowPatch.H"
#include "polyMesh.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class uniGasChapmanEnskogFreeStreamInflowFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class uniGasChapmanEnskogFreeStreamInflowFieldPatch
:
    public uniGasGeneralBoundary
{
    // Private data

        dictionary propsDict_;

        // List of (one value) mass densities corresponding to the type Ids
        List<scalarField> numberDensities_;

        // Translational temperature
        scalarField inletTransT_;

        // Rotational temperature
        scalarField inletRotT_;

        // Vibrational temperature
        scalarField inletVibT_;

        // Electronic temperature
        scalarField inletElecT_;

        // Inlet velocity
        vectorField inletVelocities_;

        // Inlet heat flux
        vectorField inletHeatFluxes_;

        // Inlet stress
        tensorField inletStresses_;

        //- List of boundary densities (one entry for each species)
        List<autoPtr<volScalarField>> boundaryNumberDensity_;

        //- Boundary temperatures
        volScalarField boundaryTransT_;

        volScalarField boundaryRotT_;

        volScalarField boundaryVibT_;

        volScalarField boundaryElecT_;

        //- Boundary velocity
        volVectorField boundaryU_;

        //- Boundary heat fluxes
        volVectorField boundaryHeatFlux_;

        //- Boundary stresses
        volTensorField boundaryStress_;

public:

    //- Runtime type information
    TypeName("uniGasChapmanEnskogFreeStreamInflowFieldPatch");

    //- Construct from components
    uniGasChapmanEnskogFreeStreamInflowFieldPatch
    (
        const polyMesh& mesh,
        uniGasCloud& cloud,
        const dictionary& dict
    );

    // Destructor
    virtual ~uniGasChapmanEnskogFreeStreamInflowFieldPatch() = default;


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
