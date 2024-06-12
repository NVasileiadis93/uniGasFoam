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
    Foam::uniGasCLLWallFieldPatch

Description
    Wall interaction based on the Cercignani-Lampis-Lord model

    Use in the boundariesDict file.

    FIELD VERSION: require boundaryU and boundaryT fields to be defined in
    your time directory.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    patchBoundaryProperties
    {
        patch       wall;
    }

    boundaryModel   uniGasCLLWallFieldPatch;

    uniGasCLLWallFieldPatchProperties
    {
        normalAccommCoeff       0.9;
        tangentialAccommCoeff   0.9;
        rotEnergyAccommCoeff    0.9;
    }
    \endverbatim

SourceFiles
    uniGasCLLWallFieldPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uniGasCLLWallFieldPatch_H
#define uniGasCLLWallFieldPatch_H

#include "uniGasPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class uniGasCLLWallFieldPatch Declaration
\*---------------------------------------------------------------------------*/

class uniGasCLLWallFieldPatch
:
    public uniGasPatchBoundary
{
    // Private data

        dictionary propsDict_;

        void setProperties();

        scalar normalAccommCoeff_;

        scalar tangentialAccommCoeff_;

        scalar rotEnergyAccommCoeff_;

        //- Boundary temperature
        volScalarField boundaryT_;

        //- Boundary velocity
        volVectorField boundaryU_;


public:

    //- Runtime type information
    TypeName("uniGasCLLWallFieldPatch");

    //- Construct from components
    uniGasCLLWallFieldPatch
    (
        const polyMesh& mesh,
        uniGasCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uniGasCLLWallFieldPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- calculate any required properties
        virtual void calculateProperties();

        //- control the mols
        virtual void controlParticle(uniGasParcel& p, uniGasParcel::trackingData& td);

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
