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
    BGKRelaxation

Description
    Relaxes particles to local Chapman-Enskog distributions

SourceFiles
    BGKRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef BGKRelaxation_H
#define BGKRelaxation_H

#include "relaxationModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "uspCloud.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        class BGKRelaxation Declaration
\*---------------------------------------------------------------------------*/

class BGKRelaxation
:
    public relaxationModel
{
    // Private data

        label infoCounter_;
        label steps_;

        scalar theta_ = 1.0;

        volScalarField relaxationFrequency_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField translationalT_;
        volVectorField UMean_;
        

        labelList typeIds_;
        scalarField rhoNMean_;
        scalarField rhoNMeanXnParticle_;
        scalarField rhoMMean_;
        scalarField rhoMMeanXnParticle_;
        scalarField linearKEMean_;
        scalarField linearKEMeanXnParticle_;
        vectorField momentumMean_;
        vectorField momentumMeanXnParticle_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("BGKRelaxation");

    //- Construct from components
    BGKRelaxation
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~BGKRelaxation() = default;


    // Member Functions

        //- Calculate macroscopic quantities
        virtual void calculateProperties();

        //- Relax particles
        virtual void relax();

        //- Sample linear velocity from CE distribution
        vector samplePostRelaxationVelocity
        (   
            scalar m,
            scalar T,
            vector U
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

