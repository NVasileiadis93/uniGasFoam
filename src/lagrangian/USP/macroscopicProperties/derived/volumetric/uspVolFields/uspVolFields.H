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
    Foam::uspVolFields

Description
    Class for measuring volumetric fields. Internal fields are measured in 
    conjunction with the cellMeasurements class. Boundary fields are measured 
    in conjunction with the boundaryMeasurements class.

    Use in fieldPropertiesDict.

Usage
    Example of field specification

    /verbatim
    field
    {
        fieldModel          	uspVolFields;
 
        timeProperties
        {
            sampleInterval              1;
            resetAtOutput               on;
 	   	    resetAtOutputUntilTime      1e-3;
        }
 
        uspVolFieldsProperties
        {
            field	                    Ar;
            typeIds                     (Ar);
            measureMeanFreePath         true;
            Tref                        1000;
            measureErrors				true;
            averagingAcrossManyRuns     true;
        }
    }
    /endverbatim

SourceFiles
    uspVolFields.C

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
        word fieldName_;
        labelList typeIds_;
        scalar Tref_;
        bool densityOnly_;
        bool measureMeanFreePath_;
        bool measureErrors_;
        bool averagingAcrossManyRuns_;

        label sampleCounter_;
        label nAvTimeSteps_;
        scalar timeAvCounter_;
        
        //- Wall unit vectors - temporary storage, should be done in boundaries/patch
        List<vectorField> n_;
        List<vectorField> t1_;
        List<vectorField> t2_;

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
        
        // Cell measurements
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

