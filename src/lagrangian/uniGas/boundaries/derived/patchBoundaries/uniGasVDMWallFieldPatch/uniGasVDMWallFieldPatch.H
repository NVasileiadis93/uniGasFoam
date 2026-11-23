/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::uniGasVDMWallFieldPatch

Description
    Wall interaction that performs a velocity dependent maxwell reflection

    Use in the boundariesDict file.

    FIELD VERSION: require boundaryU and boundaryT fields to be defined in
    your time directory.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    patchBoundaryProperties
    {
        patch               wall;
    }

    boundaryModel   uniGasVDMWallFieldPatch;

    uniGasVDMWallFieldPatchProperties
    {
        accomCoeff	    0.9;
	roughnessDegree	    0.5;
	activationEnergy    0;
	impactCoeff	    0;
    }
    \endverbatim

SourceFiles
    uniGasVDMWallFieldPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uniGasVDMWallFieldPatch_H
#define uniGasVDMWallFieldPatch_H

#include "uniGasDiffuseWallPatch.H"
#include "uniGasSpecularWallPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class uniGasVDMWallFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class uniGasVDMWallFieldPatch
:
    public uniGasPatchBoundary
{
    // Private data

        dictionary propsDict_;

	// Phi_0 in VDM paper
        scalar accomCoeff_;

        volScalarField boundaryT_;
        volVectorField boundaryU_;
	scalar temperature_;
	vector velocity_;

	// gamma in VDM paper
	scalar roughnessDegree_;

	// epsilon in VDM paper
	scalar activationEnergy_;

	// alpha in VDM orig paper
	scalar impactCoeff_;

public:

    //- Runtime type information
    TypeName("uniGasVDMWallFieldPatch");

    // Constructors

        //- Construct from components
        uniGasVDMWallFieldPatch
        (
            const polyMesh& mesh,
            uniGasCloud& cloud,
            const dictionary& dict
        );

    //- Destructor
    virtual  ~uniGasVDMWallFieldPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control the mols
        virtual void controlParticle(uniGasParcel& p, uniGasParcel::trackingData& td);

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
