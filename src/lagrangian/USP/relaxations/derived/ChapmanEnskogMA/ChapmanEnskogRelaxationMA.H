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
    ChapmanEnskogRelaxationMA

Description
    Relaxes particles to local Chapman-Enskog distributions

SourceFiles
    ChapmanEnskogRelaxationMA.C

\*---------------------------------------------------------------------------*/

#ifndef ChapmanEnskogRelaxationMA_H
#define ChapmanEnskogRelaxationMA_H

#include "relaxationModel.H"
#include "fvc.H"
#include "fvcSmooth.H"
#include "vector.H"
#include "polyMesh.H"
#include "uspCloud.H"
#include "interpolationCellPoint.H"
#include "interpolationCellPointFace.H"
#include "zeroGradientFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        class ChapmanEnskogRelaxationMA Declaration
\*---------------------------------------------------------------------------*/

class ChapmanEnskogRelaxationMA
:
    public relaxationModel
{
    // Private data


        label movingAvgCounter_;

        label movingAvgSpan_;

        label nTimeSteps_;

        label infoCounter_;

        labelList typeIds_;

        // cell measurements
        scalarField rhoNMean_;
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
        vectorField momentumMean_;
        vectorField momentumMeanXnParticle_;
        List<List<scalarField>> vibrationalETotal_;
        List<scalarField> electronicETotal_;
        List<scalarField> nParcels_;
        List<scalarField> mccSpecies_;
        List<scalarField> nGroundElectronicLevel_;
        List<scalarField> nFirstElectronicLevel_;

        // boundary measurements
        List<labelList> boundaryCells_;
        
        List<scalarField> rhoNBF_;
        List<scalarField> rhoMBF_;
        List<scalarField> linearKEBF_;
        List<scalarField> rotationalEBF_;
        List<scalarField> rotationalDofBF_;
        List<scalarField> vibTxvDofBF_;
        List<scalarField> totalvDofBF_;
        List<scalarField> speciesRhoNIntBF_;
        List<scalarField> speciesRhoNElecBF_;
        List<vectorField> momentumBF_;

        List<List<scalarField>> vibrationalEBF_;
        List<List<scalarField>> electronicEBF_;
        List<List<scalarField>> speciesRhoNBF_;
        List<List<scalarField>> mccSpeciesBF_;
        List<List<scalarField>> vibTBF_;
        List<List<scalarField>> vDofBF_;

        // vol fields
        volScalarField Prandtl_;
        volScalarField relaxFreq_;
        volScalarField viscosity_;
        volScalarField conductivity_;
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
        volVectorField gradT_;
        volTensorField gradU_;

        volScalarField translationalTAvg_;
        volVectorField UMeanAvg_;
        volVectorField heatFluxVectorAvg_;
        volTensorField stressTensorAvg_;

        surfaceScalarField phi_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("ChapmanEnskogRelaxationMA");

    //- Construct from components
    ChapmanEnskogRelaxationMA
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~ChapmanEnskogRelaxationMA() = default;


    // Member Functions

        //- Calculate macroscopic quantities
        virtual void calculateProperties();

        //- Reset macroscopic quantities
        virtual void resetProperties();

        //- Change post-collision velocities to conserve momentum and energy
        virtual void conserveMomentumAndEnergy(const label& cell);

        //- Relax particles
        virtual void relax();

        //- Sample linear velocity from CE distribution
        vector samplePostRelaxationVelocity
        (   
            const scalar& breakdownParameter,
            const scalar& mass,
            const scalar& p,
            const scalar& translationalT,
            const vector& UMean,
            const vector& heatFluxVector,
            const tensor& stressTensor
        );


        //- Sample linear velocity from CE distribution
        vector samplePostRelaxationVelocityTest
        (   
            const scalar& m,
            const scalar& Pr,
            const scalar& p,
            const scalar& T,
            const vector& U,
            const tensor& pT
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

