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
    Foam::uniGasForceSurface

Description
    Class for measuring fluxes through a surface. Fields are measured in 
    conjunction with the faceTracker class. 

    Use in fieldPropertiesDict.

Usage
    Example of field specification

    /verbatim
     field
     {
         fieldModel          	        uniGasForceSurface;
 
         timeProperties
         {
            sampleInterval              1;
         	resetAtOutput               on;
         	resetAtOutputUntilTime  	1e-3;
         }
         
         uniGasForceSurfaceProperties
         {
            field                       Ar;
         	typeIds                     (Ar);
            faceZone                    sliceFaceZone;
            fluxDirection               (1 0 0);
            averagingAcrossManyRuns     true;
         }
     } 
    /endverbatim

SourceFiles
    uniGasForceSurface.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasForceSurface_H
#define uniGasForceSurface_H

#include "uniGasField.H"
#include "vector.H"
#include "uniGasCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uniGasForceSurface Declaration
\*---------------------------------------------------------------------------*/

class uniGasForceSurface
:
	public uniGasField
{
    // Private data

        //- Properties dictionary
        dictionary propsDict_;

        word fieldName_;

        labelList typeIds_;

        word patchName_;  

        bool averagingAcrossManyRuns_;

        label sampleCounter_;
        
        label patchId;

        scalar timeAvCounter_;

        label timeIndex_;

        vector force_;

        vectorField forcePatch_;

    // Private Member Functions

        void readIn();

        void writeOut();

 public:

    //- Runtime type information
    TypeName("uniGasForceSurface");

    // Constructors

        //- Construct from components
        uniGasForceSurface
        (
            const Time& t,
            const polyMesh& mesh,
            uniGasCloud& molCloud,
            const dictionary& dict
        );

    //- Destructor
    virtual ~uniGasForceSurface() = default;


    // Member Functions

        //- Create initial field
        virtual void createField();

        //- Calculate field
        virtual void calculateField();


        //- Write field
        virtual void writeField();

        virtual void updateProperties(const dictionary& dict);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
