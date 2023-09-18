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
    Foam::uspCLLWallPatch

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
        patchName                           wall;
    }

    boundaryModel   uspCLLWallPatch;

    uspCLLWallPatchProperties
    {
        temperature                                 300;
        velocity                                    (0 0 0);
        normalAccommodationCoefficient              0.9;
        tangentialAccommodationCoefficient          0.9;
        rotationalEnergyAccommodationCoefficient    0.9;
        vibrationalEnergyAccommodationCoefficient   0.9;
    }
    \endverbatim

SourceFiles
    uspCLLWallPatch.C

\*---------------------------------------------------------------------------*/

#ifndef uspCLLWallPatch_H
#define uspCLLWallPatch_H

#include "uspPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class uspCLLWallPatch Declaration
\*---------------------------------------------------------------------------*/

class uspCLLWallPatch
:
    public uspPatchBoundary
{
    // Private data

        dictionary propsDict_;

        scalar normalAccommodationCoefficient_;

        scalar tangentialAccommodationCoefficient_;

        scalar rotationalEnergyAccommodationCoefficient_;

        scalar vibrationalEnergyAccommodationCoefficient_;

        void setProperties();


public:


    //- Runtime type information
    TypeName("uspCLLWallPatch");

    //- Construct from components
    uspCLLWallPatch
    (
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspCLLWallPatch() = default;


    // Member Functions

        virtual void initialConfiguration();

        //- Calculate any required properties
        virtual void calculateProperties();

        //- Control the mols
        virtual void controlParticle
        (
            uspParcel& p,
            uspParcel::trackingData& td
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
