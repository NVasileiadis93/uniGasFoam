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
    uspVolFields

Description
    Creates geometric fields for temperature, density, velocity, pressure.

    Measures overall temperature, including vibrational temperature, for a
    single specie gas or a gas mixture and writes the results to a volume
    scalar field for further post-processing.

    Translational, rotatational and vibrational temperature fields will also be
    written automatically.

    Boundary fields are measured in conjunction with the boundaryMeasurements
    class and are also written.

SourceFiles
    uspVolFieldsI.H
    uspVolFields.C
    uspVolFieldsIO.C

\*---------------------------------------------------------------------------*/

#ifndef uspVolFields_H
#define uspVolFields_H

#include "uspField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        class uspVolFields Declaration
\*---------------------------------------------------------------------------*/

class uspVolFields
:
    public uspField
{
    // Private data

        dictionary propsDict_;
        
        //- Wall unit vectors - temporary storage, should be done in
        //  boundaries/patch
        List<vectorField> n_;
        List<vectorField> t1_;
        List<vectorField> t2_;

        bool resetFieldsAtOutput_;
        scalar resetFieldsAtOutputUntilTime_;
        label sampleInterval_;
        label sampleCounter_;
        scalar mfpReferenceTemperature_;
        word fieldName_;
        volScalarField uspRhoN_;
        volScalarField uspRhoNMean_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField translationalT_;
        volScalarField rotationalT_;
        volScalarField vibrationalT_;
        volScalarField electronicT_;
        volScalarField overallT_;
        volScalarField q_;
        volScalarField tau_;
        volScalarField MFP_; // Mean free path
        volScalarField dxMFP_;
        volScalarField MCR_; // Mean collision rate
        volScalarField MCT_;
        volScalarField dtMCT_;
        volScalarField Ma_;
        volScalarField densityError_;
        volScalarField velocityError_;
        volScalarField temperatureError_;
        volScalarField pressureError_;
        volVectorField UMean_;
        volVectorField fD_;
        volVectorField heatFluxVector_;
        volTensorField pressureTensor_;
        volTensorField shearStressTensor_;
        label nTimeSteps_;
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
        scalarField nColls_;
        vectorField momentumMean_;
        vectorField momentumMeanXnParticle_;
        List<labelList> boundaryCells_;
        List<List<scalarField>> vibrationalETotal_;
        List<scalarField> electronicETotal_;
        List<scalarField> nParcels_;
        List<scalarField> nParcelsXnParticle_;
        List<scalarField> mccSpecies_;
        List<scalarField> nGroundElectronicLevel_;
        List<scalarField> nFirstElectronicLevel_;
        List<scalarField> mfp_; // mean free path
        List<scalarField> mcr_; // mean collision rate


        // boundary measurements

        List<scalarField> rhoNBF_;
        List<scalarField> rhoMBF_;
        List<scalarField> linearKEBF_;
        List<scalarField> rotationalEBF_;
        List<scalarField> rotationalDofBF_;
        List<scalarField> qBF_;
        List<scalarField> totalvDofBF_;
        List<scalarField> speciesRhoNIntBF_;
        List<scalarField> speciesRhoNElecBF_;
        List<vectorField> momentumBF_;
        List<vectorField> fDBF_;

        List<List<scalarField>> vibrationalEBF_;
        List<List<scalarField>> electronicEBF_;
        List<List<scalarField>> speciesRhoNBF_;
        List<List<scalarField>> mccSpeciesBF_;

        bool averagingAcrossManyRuns_;
        bool measureMeanFreePath_;
        bool measureErrors_;
        bool densityOnly_;
        bool measureHeatFluxShearStress_;

    // Private Member Functions

        void readIn();

        void writeOut();

        void calculateWallUnitVectors();


public:

    //- Runtime type information
    TypeName("uspVolFields");

    //- Construct from components
    uspVolFields
    (
        const Time& t,
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspVolFields() = default;


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

