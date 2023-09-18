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
    ChapmanEnskogRelaxation

Description
    Relaxes particles to local Chapman-Enskog distributions

SourceFiles
    ChapmanEnskogRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef ChapmanEnskogRelaxation_H
#define ChapmanEnskogRelaxation_H

#include "relaxationModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "uspCloud.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        class ChapmanEnskogRelaxation Declaration
\*---------------------------------------------------------------------------*/

class ChapmanEnskogRelaxation
:
    public relaxationModel
{
    // Private data

        label infoCounter_;
        scalar Pr_=1.0;
        
        volScalarField relaxationFrequency_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField translationalT_;
        volScalarField rotationalT_;
        volScalarField vibrationalT_;
        volScalarField electronicT_;
        volScalarField overallT_;
        volVectorField UMean_;
        volVectorField heatFluxVector_;
        volTensorField pressureTensor_;
        volTensorField stressTensor_;

        labelList typeIds_;
        scalarField rhoNMean_;
        scalarField rhoNInstantaneous_;
        scalarField rhoNMeanXnParticle_;
        scalarField rhoNMeanInt_;
        scalarField molsElec_;
        scalarField rhoMMean_;
        scalarField rhoMMeanXnParticle_;
        scalarField linearKEMean_;
        scalarField linearKEMeanXnParticle_;
        scalarField rotationalEMean_;
        scalarField rotationalDofMean_;
        scalarField muu_;
        scalarField muv_;
        scalarField muw_;
        scalarField mvv_;
        scalarField mvw_;
        scalarField mww_;
        scalarField mcc_;
        scalarField mccu_;
        scalarField mccv_;
        scalarField mccw_;
        scalarField eu_;
        scalarField ev_;
        scalarField ew_;
        scalarField e_;
        scalarField totalvDof_;
        vectorField momentumMean_;
        vectorField momentumMeanXnParticle_;
        List<List<scalarField>> vibrationalETotal_;
        List<scalarField> electronicETotal_;
        List<scalarField> nParcels_;
        List<scalarField> nParcelsXnParticle_;
        List<scalarField> mccSpecies_;
        List<scalarField> vibT_;
        List<scalarField> vDof_;
        List<scalarField> nGroundElectronicLevel_;
        List<scalarField> nFirstElectronicLevel_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("ChapmanEnskogRelaxation");

    //- Construct from components
    ChapmanEnskogRelaxation
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~ChapmanEnskogRelaxation() = default;


    // Member Functions

        //- Calculate macroscopic quantities
        virtual void calculateProperties();

        //- Relax particles
        virtual void relax();

        //- Sample linear velocity from CE distribution
        vector samplePostRelaxationVelocity
        (   
            const scalar& B,
            const scalar& m,
            const scalar& p,
            const scalar& T,
            const vector& U,
            const vector& q,
            const tensor& s
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

