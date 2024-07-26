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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::uniGasReflectiveParticleMembranePatch

Description
    This class provides functions to implement a cyclic reflective membrane.

    Use in the boundariesDict file.

Usage
    Example of the boundary condition specification:
    
    \verbatim
    generalBoundaryProperties
    {
        patch       membrane;
    }

    boundaryModel   uniGasReflectiveParticleMembranePatch;

    uniGasReflectiveParticleMembranePatchProperties
    {
        reflectionProbability       0.5;
        temperature                 300;
        velocity                    (0 0 0);
    }
    \endverbatim

SourceFiles
    uniGasReflectiveParticleMembranePatch.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasReflectiveParticleMembranePatch_H
#define uniGasReflectiveParticleMembranePatch_H

#include "uniGasCyclicBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class uniGasReflectiveParticleMembranePatch Declaration
\*---------------------------------------------------------------------------*/

class uniGasReflectiveParticleMembranePatch
:
	public uniGasCyclicBoundary
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
    TypeName("uniGasReflectiveParticleMembranePatch");

    // Constructors

        //- Construct from components
        uniGasReflectiveParticleMembranePatch
		(
            const polyMesh& mesh,
            uniGasCloud& cloud,
            const dictionary& dict
		);


    //- Destructor
    virtual ~uniGasReflectiveParticleMembranePatch() = default;


    // Member Functions

        void calculateProperties();

        //- Create an initial configuration
        void initialConfiguration();

        //- Control the molecules
        void controlMol
        (
            uniGasParcel& mol,
            uniGasParcel::trackingData& td
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
