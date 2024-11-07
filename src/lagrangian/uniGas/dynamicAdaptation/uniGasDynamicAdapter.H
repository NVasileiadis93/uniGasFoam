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
    Foam::uniGasDynamicAdapter

Description
    Class providing functions for dynamically adapting cloud quantities

SourceFiles
    uniGasDynamicAdapter.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasDynamicAdapter_H
#define uniGasDynamicAdapter_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "fvc.H"
#include "fvcSmooth.H"
#include "volFields.H"
#include "vectorField.H"
#include "tensorField.H"
#include "Random.H"
#include "uniGasParcel.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uniGasDynamicAdapter Declaration
\*---------------------------------------------------------------------------*/

class uniGasDynamicAdapter
{

    // Private data

        const dictionary& dict_;

        Time& time_;

        fvMesh& mesh_;

        uniGasCloud& cloud_;

        const label minSubCellLevels_;

        const label maxSubCellLevels_;

        bool timeStepAdaptation_;

        bool subCellAdaptation_;

        bool cellWeightAdaptation_;

        label adaptationInterval_;

        scalar maxTimeStepMCTRatio_;

        scalar maxCourantNumber_;

        scalar maxSubCellSizeMFPRatio_;

        label smoothingPasses_;

        scalar theta_;

        label timeSteps_;

        scalar timeAvCounter_;

        labelList typeIds_;        

        // cell measurements
        scalarField rhoNMean_;
        scalarField rhoNMeanXnParticle_;
        scalarField rhoMMeanXnParticle_;
        scalarField linearKEMeanXnParticle_;
        vectorField momentumMeanXnParticle_;

        List<scalarField> nParcelsXnParticle_;
        
        // fields required for dynamic adaptation
        volScalarField rhoN_;
        volScalarField translationalT_;
        volVectorField UMean_;

        volScalarField cellWeightFactor_;
        volScalarField timeStepMCTRatio_;
        volVectorField courantNumber_;
        volVectorField cellSizeMFPRatio_;
        volVectorField prevCellSizeMFPRatio_;

    // Private Member Functions

public:

    //- Construct from components
    uniGasDynamicAdapter
    (
        const dictionary& dict,
        Time& time,
        fvMesh& mesh,
        uniGasCloud& cloud
    );

    //- Destructor
    virtual ~uniGasDynamicAdapter() = default;

    // Member Functions
    virtual void setInitialConfiguration
    (
        const Field<scalar>& speciesRhoN,
        const scalar& transT,
        const vector& U
    );

    virtual void setInitialConfiguration
    (
        const cellZone& zone,
        const Field<scalar>& speciesRhoN,
        const scalar& transT,
        const vector& U
    );

    virtual void setInitialConfiguration
    (
        const List<autoPtr<volScalarField>>& speciesRhoNPtr,
        const volScalarField& transT,
        const volVectorField& U
    );
    
    void calculateAdaptationQuantities
    (
        const label& cell,
        const scalar& rhoN,
        const scalar& transT,
        const vector& U,
        const scalarList& speciesRhoN,
        scalar& timeStepMCTRatio,
        vector& courantNumber,
        vector& cellSizeMFPRatio
    );

    void calculateTimeStep();

    vector calculateSubCellLevels
    (
        const label& cell,
        const vector& cellSizeMFPRatio
    );

    scalar calculateCellWeightFactor
    (
        const label& cell,
        const scalar& rhoN
    );

    virtual void smoothCellWeightFactor
    (
        volScalarField& rhoN,
        volScalarField& cellWeightFactor
    );

    void adapt();

    // Acess
    inline const label& minSubCellLevels() const;

    inline const label& maxSubCellLevels() const;

    inline const bool& timeStepAdaptation() const;

    inline const bool& cellWeightAdaptation() const;

    inline const bool& subCellAdaptation() const;

    inline const scalar& maxSubCellSizeMFPRatio() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uniGasDynamicAdapterI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

