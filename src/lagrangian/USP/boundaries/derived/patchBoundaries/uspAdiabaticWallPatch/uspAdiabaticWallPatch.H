/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    uspAdiabaticWallPatch

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    SINGLE VALUE: temperature and velocity required in boundariesDict


SourceFiles
    uspAdiabaticWallPatch.C


\*---------------------------------------------------------------------------*/

#ifndef uspAdiabaticWallPatch_H
#define uspAdiabaticWallPatch_H

#include "uspPatchBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "uspCloud.H"
#include "distribution.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uspAdiabaticWallPatch Declaration
\*---------------------------------------------------------------------------*/

class uspAdiabaticWallPatch
:
	public uspPatchBoundary
{

private:

    // Private data
        dictionary propsDict_;
        vector incidentMomentum_;
        vector gasVelocity_;
        scalar totalVolume_;
        scalar totalMass_;
        scalar incidentMass_;
        scalar wallTemperature_;
        scalar mols_;
        scalar totalMols_;
        scalar cv_;
        scalar mcc_;

        void setProperties();

public:

    //- Runtime type information
    TypeName("uspAdiabaticWallPatch");

    // Constructors

        //- Construct from components
        uspAdiabaticWallPatch
        (
            Time& t,
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        );

    // Destructor

        ~uspAdiabaticWallPatch();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        void controlParticle(uspParcel& p, uspParcel::trackingData& td);

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );
        
        void writeField();

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);



        // Access





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
