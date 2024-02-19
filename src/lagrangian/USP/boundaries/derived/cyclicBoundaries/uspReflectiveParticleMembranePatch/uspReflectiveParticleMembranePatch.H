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
    Foam::uspReflectiveParticleMembranePatch

Description
    RPM technique Li et al.

SourceFiles
    uspReflectiveParticleMembranePatch.C

\*---------------------------------------------------------------------------*/

#ifndef uspReflectiveParticleMembranePatch_H
#define uspReflectiveParticleMembranePatch_H

#include "uspCyclicBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class uspReflectiveParticleMembranePatch Declaration
\*---------------------------------------------------------------------------*/

class uspReflectiveParticleMembranePatch
:
	public uspCyclicBoundary
{
    // Private data

        dictionary propsDict_;

        scalar p_;

        scalar temperature_;

        vector velocity_;

        label nReflections_;
        
        label nRejections_;

public:

    //- Runtime type information
    TypeName("uspReflectiveParticleMembranePatch");

    // Constructors

        //- Construct from components
        uspReflectiveParticleMembranePatch
		(
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
		);


    //- Destructor
    virtual ~uspReflectiveParticleMembranePatch() = default;


    // Member Functions

        void calculateProperties();

        //- Create an initial configuration
        void initialConfiguration();

        //- Control the molecules
        void controlMol
        (
            uspParcel& mol,
            uspParcel::trackingData& td
        );

        //- Output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- Update properties from a modified dictionary
        void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
