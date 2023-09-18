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
    uspBinsMethod

Description


SourceFiles
    uspBinsMethodI.H
    uspBinsMethod.C
    uspBinsMethodIO.C

\*---------------------------------------------------------------------------*/

#ifndef uspBinsMethod_H
#define uspBinsMethod_H

#include "uspField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "uspCloud.H"
#include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uspBinsMethod Declaration
\*---------------------------------------------------------------------------*/

class uspBinsMethod
:
	public uspField
{
    // Private data

        //- subDictionary containing the properties
        dictionary propsDict_;

        autoPtr<binModel> binModel_;

        word regionName_;
        label regionId_;
        word fieldName_;

        labelList typeIds_;

        label sampleInterval_;
        label sampleCounter_;

        scalar averagingCounter_;

        scalarField volumeOfCellsInBin_;
        scalarField mols_;
        scalarField uspMols_;
        scalarField molsInt_;
        scalarField molsElec_;
        scalarField mass_;
        scalarField mcc_;
        vectorField mom_;
        vectorField UCollected_;
        scalarField rotationalEMean_;
        scalarField rotationalDofMean_;

        // pressure measurements

        scalarField muu_;
        scalarField muv_;
        scalarField muw_;
        scalarField mvv_;
        scalarField mvw_;
        scalarField mww_;
        scalarField scalarPressure_;
        scalarField mccu_;
        scalarField mccv_;
        scalarField mccw_;

        scalarField eu_;
        scalarField ev_;
        scalarField ew_;
        scalarField e_;

        List<scalarField> speciesMols_;
        List<scalarField> mccSpecies_;
        List<scalarField> electronicETotal_;
        List<scalarField> nParticlesGroundElectronicState_;
        List<scalarField> nParticlesFirstElectronicState_;
        List<scalarField> vDof_;
        List<scalarField> mfp_;

        scalarField binVolume_;
        scalarField N_;
        scalarField M_;
        scalarField rhoN_;
        scalarField rhoM_;
        scalarField averageMass_;
        vectorField UMean_;
        vectorField UCAM_;
        scalarField scalarUMean_;
        scalarField translationalTemperature_;
        scalarField rotationalTemperature_;
        scalarField vibrationalTemperature_;
        scalarField electronicTemperature_;
        scalarField overallTemperature_;
        tensorField pField_;
        tensorField tauField_;
        vectorField qField_;
        vectorField qInternalField_;
        vectorField qTranslationalField_;
        scalarField meanFreePath_;
        scalarField Ma_;

        List<List<scalarField>> vibrationalETotal_;

        List<bool> outputField_;

        bool averagingAcrossManyRuns_;
        bool permeabilityMeasurements_;
        bool rotationalMeasurements_;
        scalar wallVelocity_;
        scalar wallRadius_;


    // Private Member Functions

        void readIn();

        void writeOut();


public:

    //- Runtime type information
    TypeName("uspBinsMethod");

    // Constructors

        //- Construct from components
        uspBinsMethod
        (
            const Time& t,
            const polyMesh& mesh,
            uspCloud& cloud,
            const dictionary& dict
        );

    // Destructor
    virtual ~uspBinsMethod() = default;


    // Member Functions

        //- create initial field
        virtual void createField();

        //- calculate field
        virtual void calculateField();

        //- return field
//         const volScalarField& densityField() const;
//         const propertyField& fields() const;
        //- write field
        virtual void writeField();

        virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
