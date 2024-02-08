/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Foam::uspDynamicAdapter

Description
    Basic/abstract class of a controller

SourceFiles
    uspDynamicAdapter.C

\*---------------------------------------------------------------------------*/

#ifndef uspDynamicAdapter_H
#define uspDynamicAdapter_H

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
#include "uspParcel.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uspDynamicAdapter Declaration
\*---------------------------------------------------------------------------*/

class uspDynamicAdapter
{

    // Private data

        label minSubcellLevels_ = 2;

        label maxSubcellLevels_ = 5;

        label smoothingPasses_;

        scalar theta_;

        const dictionary& dict_;

        const fvMesh& mesh_;

        uspCloud& cloud_;

        Random& rndGen_;

        label timeSteps_;

        label nAvTimeSteps_;

        bool adaptiveSubcells_;

        bool adaptiveCellWeights_;

        label timeInterval_;

        scalar maxSubcellSizeMFPRatio_;

        scalar Tref_;

        labelList typeIds_;        

        // cell measurements
        scalarField rhoNMean_;
        scalarField rhoNMeanXnParticle_;
        scalarField rhoMMeanXnParticle_;
        scalarField linearKEMeanXnParticle_;
        vectorField momentumMeanXnParticle_;
        scalarField MFP_;
        scalarField MCT_;

        List<scalarField> nParcels_;
        List<scalarField> nParcelsXnParticle_;
        List<scalarField> mfp_;
        List<scalarField> mct_;
        
        // fields required for dynamic adaptation
        volScalarField rhoN_;
        volScalarField translationalT_;
        volVectorField UMean_;
        volVectorField cellSizeMFPRatio_;
        volScalarField timeStepMCTRatio_;
        volScalarField cellWeightFactor_;

    // Private Member Functions

public:

    //- Construct from components
    uspDynamicAdapter
    (
        const dictionary& dict,
        const fvMesh& mesh,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~uspDynamicAdapter() = default;

    // Member Functions

    void update();
 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "uspDynamicLoadBalancingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

