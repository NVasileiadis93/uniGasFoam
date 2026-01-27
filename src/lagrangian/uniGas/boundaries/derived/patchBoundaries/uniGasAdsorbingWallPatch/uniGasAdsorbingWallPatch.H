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
    Foam::uniGasAdsorbingWallPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature.

    Use in the boundariesDict file.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    patchBoundaryProperties
    {
        patch       wall;
    }

    boundaryModel   uniGasAdsorbingWallPatch;

    uniGasAdsorbingWallPatchProperties
    {
        velocity 	    		(0 0 0);
        temperature     		273;
		adsorbedIds				(Ar H2);
		adsorptionProbabilities {Ar 0.5; H2 0.8;};
    }
    \endverbatim


SourceFiles
    uniGasAdsorbingWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uniGasAdsorbingWallPatch_H
#define uniGasAdsorbingWallPatch_H

#include "uniGasPatchBoundary.H"
#include "uniGasCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class uniGasAdsorbingWallPatch Declaration
\*---------------------------------------------------------------------------*/

class uniGasAdsorbingWallPatch
:
    public uniGasPatchBoundary
{
    // Private data

        dictionary propsDict_;

         //- typeIds to be adsorbed
        labelList typeIds_;
        
        //- probability of adsorbtion for each species in typeIds_
        scalarList adsorptionProbs_;

        //- wall temperature
        scalar temperature_;

        //- wall velocity
        vector velocity_;

        void setProperties();

public:

    //- Runtime type information
    TypeName("uniGasAdsorbingWallPatch");

    //- Construct from components
    uniGasAdsorbingWallPatch
    (
        const polyMesh& mesh,
        uniGasCloud& cloud,
        const dictionary& dict
    );

    // Destructor
    virtual ~uniGasAdsorbingWallPatch() = default;


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
