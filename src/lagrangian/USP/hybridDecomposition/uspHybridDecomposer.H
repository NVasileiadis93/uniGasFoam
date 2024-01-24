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
    Foam::uspHybridDecomposer

Description
    Basic/abstract class of a controller

SourceFiles
    uspHybridDecomposer.C

\*---------------------------------------------------------------------------*/

#ifndef uspHybridDecomposer_H
#define uspHybridDecomposer_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "fvc.H"
#include "volFields.H"
#include "vectorField.H"
#include "tensorField.H"
#include "SortableList.H"
#include "Random.H"
#include "uspParcel.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uspHybridDecomposer Declaration
\*---------------------------------------------------------------------------*/

class uspHybridDecomposer
{

    // Private data

        const dictionary& dict_;

        const fvMesh& mesh_;

        uspCloud& cloud_;

        Random& rndGen_;

        bool firstDecomp;

        label timeSteps_;

        label nAvTimeSteps_;

        word breakdownCriterion_;

        label decomposeInterval_;

        scalar breakdownMax_;

        scalar theta_;

        label smoothingPasses_;

        label refinementPasses_;

        scalar Tref_;

        labelList typeIds_;        

        // cell measurements
        scalarField rhoNMean_;
        scalarField rhoNMeanXnParticle_;
        scalarField rhoMMean_;
        scalarField rhoMMeanXnParticle_;
        scalarField linearKEMean_;
        scalarField linearKEMeanXnParticle_;
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
        List<scalarField> nParcels_;
        List<scalarField> nParcelsXnParticle_;
        List<scalarField> mccSpecies_;

        // boundary measurements
        List<labelList> boundaryCells_;
        
        List<scalarField> rhoNBF_;
        List<scalarField> rhoMBF_;
        List<scalarField> linearKEBF_;
        List<vectorField> momentumBF_;
        List<List<scalarField>> speciesRhoNBF_;
        List<List<scalarField>> mccSpeciesBF_;

        // fields required for decomposition
        volScalarField B_;
        volScalarField KnGLL_;
        volScalarField KnRho_;
        volScalarField CLB_;
        volVectorField CLBQ_;
        volTensorField CLBS_;
        volScalarField ASq_;
        volScalarField breakdownParameter_;
        volScalarField KnT_;
        volScalarField KnU_;
        volScalarField mu_;
        volScalarField kappa_;
        volScalarField prevRhoN_;
        volScalarField prevRhoM_;
        volScalarField prevP_;
        volScalarField prevTranslationalT_;
        volVectorField prevUMean_;
        volVectorField prevHeatFluxVector_;
        volTensorField prevPressureTensor_;
        volTensorField prevShearStressTensor_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField translationalT_;
        volVectorField UMean_;
        volVectorField heatFluxVector_;
        volTensorField pressureTensor_;
        volTensorField shearStressTensor_;
        volVectorField heatFluxVectorNS_;
        volTensorField shearStressTensorNS_;
    // Private Member Functions

public:

    //- Construct from components
    uspHybridDecomposer
    (
        const dictionary& dict,
        const fvMesh& mesh_,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~uspHybridDecomposer() = default;

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

