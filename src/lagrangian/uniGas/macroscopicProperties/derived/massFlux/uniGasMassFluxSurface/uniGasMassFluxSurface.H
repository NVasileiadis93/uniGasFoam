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
    Foam::uniGasMassFluxSurface

Description
    Class for measuring fluxes through a surface. Fields are measured in 
    conjunction with the faceTracker class. 

    Use in fieldPropertiesDict.

Usage
    Example of field specification

    /verbatim
     field
     {
         fieldModel          	        uniGasMassFluxSurface;
 
         timeProperties
         {
            sampleInterval              1;
         	resetAtOutput               on;
         	resetAtOutputUntilTime  	1e-3;
         }
         
         uniGasMassFluxSurfaceProperties
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
    uniGasMassFluxSurface.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasMassFluxSurface_H
#define uniGasMassFluxSurface_H

#include "uniGasField.H"
#include "vector.H"
#include "uniGasCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uniGasMassFluxSurface Declaration
\*---------------------------------------------------------------------------*/

class uniGasMassFluxSurface
:
	public uniGasField
{
    // Private data

        //- Properties dictionary
        dictionary propsDict_;

        word fieldName_;

        labelList typeIds_;

        word faceZoneName_;  

        vector fluxDirection_;

        bool averagingAcrossManyRuns_;

        label sampleCounter_;
        
        label regionId_;

        scalar zoneSurfaceArea_;

        scalar molsZone_;

        scalar massZone_;

        scalar momentumZone_;

        scalar energyZone_;

        scalar timeAvCounter_;

        label timeIndex_;

        scalarField molFluxZone_;
        scalarField massFluxZone_;
        scalarField momentumFluxZone_;
        scalarField energyFluxZone_;

    // Private Member Functions

        void readIn();

        void writeOut();


 public:

    //- Runtime type information
    TypeName("uniGasMassFluxSurface");

    // Constructors

        //- Construct from components
        uniGasMassFluxSurface
        (
            const Time& t,
            const polyMesh& mesh,
            uniGasCloud& molCloud,
            const dictionary& dict
        );

    //- Destructor
    virtual ~uniGasMassFluxSurface() = default;


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
