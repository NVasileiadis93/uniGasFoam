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
    Foam::uniGasCLLWallPatch

Description
    Cercignani-Lampis-Lord surface interaction model - currently deals with
    the translation and rotational modes. No vibrational energy has been dealt
    with yet and it does not include diffuse scattering with incomplete
    energy accommodation.

    Use in the boundariesDict file.

    References:
        C. Cercignani and M. Lampis, " Kinetic models for gas-surface interactions",
        Transport Theory and Statistical Physics, 1(2), 101-114, 1971.

        R.G. Lord, "Some extensions to the Cercignani–Lampis gas–surface scattering kernel",
        Physics of Fluids A, 3(4), 706-710, 1991.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    patchBoundaryProperties
    {
        patch       wall;
    }

    boundaryModel   uniGasCLLWallPatch;

    uniGasCLLWallPatchProperties
    {
        normalAccommCoeff       0.9;
        tangentialAccommCoeff   0.9;
        rotEnergyAccommCoeff    0.9;
        temperature             300;
        velocity                (0 0 0);
    }
    \endverbatim

SourceFiles
    uniGasCLLWallPatch.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasCLLWallPatch_H
#define uniGasCLLWallPatch_H

#include "uniGasPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class uniGasCLLWallPatch Declaration
\*---------------------------------------------------------------------------*/

class uniGasCLLWallPatch
:
    public uniGasPatchBoundary
{
    // Private data

        dictionary propsDict_;

        scalar normalAccommCoeff_;

        scalar tangentialAccommCoeff_;

        scalar rotEnergyAccommCoeff_;

        scalar vibEnergyAccommCoeff_;

        scalar temperature_;

        vector velocity_;

public:


    //- Runtime type information
    TypeName("uniGasCLLWallPatch");

    //- Construct from components
    uniGasCLLWallPatch
    (
        const polyMesh& mesh,
        uniGasCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uniGasCLLWallPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control the mols
        virtual void controlParticle
        (
            uniGasParcel& p,
            uniGasParcel::trackingData& td
        );

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
