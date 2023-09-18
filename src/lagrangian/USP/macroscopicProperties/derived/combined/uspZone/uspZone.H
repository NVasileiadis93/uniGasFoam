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
    uspZone

Description


SourceFiles
    uspZoneI.H
    uspZone.C
    uspZoneIO.C

\*---------------------------------------------------------------------------*/

#ifndef uspZone_H
#define uspZone_H

#include "uspField.H"
#include "vector.H"
#include "polyMesh.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uspZone Declaration
\*---------------------------------------------------------------------------*/

class uspZone
:
	public uspField
{
    // Private data

        //- Properties dictionary
        dictionary propsDict_;

        word regionName_;
        label regionId_;
        label sampleInterval_;
        label sampleCounter_;
        word fieldName_;
        scalar totalVolume_;

        labelList typeIds_;
        label timeIndex_;
        scalar averagingCounter_;

        scalar mols_;
        scalar uspMols_;
        scalar molsInt_;
        scalar molsElec_;
        scalar mass_;
        scalar mcc_;
        vector mom_;
        vector UCollected_;
        scalar rotationalEMean_;
        scalar rotationalDofMean_;

        symmTensor muu_;

        vector mccu_;
        scalar nColls_;

        vector eu_;
        scalar e_;

        label stepIndex_;
        scalarField speciesMols_;
        scalarField mccSpecies_;
        scalarField electronicETotal_;
        scalarField nParticlesGroundElectronicState_;
        scalarField nParticlesFirstElectronicState_;
        scalarField vDof_;
        scalarField mfp_;
        scalarField mcr_; // mean collision rate
        label nSteps_;

        scalarField N_;
        scalarField rhoN_;
        scalarField rhoM_;
        vectorField UMean_;
        vectorField UCAM_;
        scalarField translationalTemperature_;
        scalarField rotationalTemperature_;
        scalarField vibrationalTemperature_;
        scalarField electronicTemperature_;
        scalarField overallTemperature_;
        scalarField scalarPressure_;
        symmTensorField pField_;
        symmTensorField tauField_;
        vectorField qField_;
        vectorField qInternalField_;
        vectorField qTranslationalField_;
        vectorField vibrationalModeTemperatures_;
        scalarField meanFreePath_;
        scalarField meanCollisionRate_;
        scalarField meanCollisionTime_;
        scalarField meanCollisionTimeTimeStepRatio_;
        scalarField measuredCollisionRate_;
        scalarField Ma_;

        List<scalarField> vibrationalETotal_;

        List<bool> outputField_;

        bool instantaneous_;

        bool averagingAcrossManyRuns_;

        bool measureCollisionRateOnly_;


    // Private Member Functions

        void read();

        void write();


public:

    //- Runtime type information
    TypeName("uspZone");

    //- Construct from components
    uspZone
    (
        const Time& t,
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspZone() = default;


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
