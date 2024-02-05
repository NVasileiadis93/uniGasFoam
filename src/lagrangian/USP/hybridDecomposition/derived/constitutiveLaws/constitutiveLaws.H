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
    Foam::constitutiveLaws

Description
    Hybrid decomposition class based on constitutive law breakdown criterion

SourceFiles
    constitutiveLaws.C

\*---------------------------------------------------------------------------*/

#ifndef constitutiveLaws_H
#define constitutiveLaws_H

#include "uspHybridDecomposition.H"
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
#include "uspCloud.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class constitutiveLaws Declaration
\*---------------------------------------------------------------------------*/

class constitutiveLaws
:
    public uspHybridDecomposition
{

    // Private data

        label decomposeInterval_;
        scalar breakdownMax_;
        bool timeAverage_;
        scalar theta_;
        label smoothingPasses_;
        scalar Tref_;

        bool firstDecomp_;
        label timeSteps_;
        label nAvTimeSteps_;

        label refinementPasses_ = 10;
        labelList typeIds_;        

        // cell measurements
        scalarField rhoNMean_;
        scalarField rhoMMean_;
        scalarField linearKEMean_;
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
        scalarField rhoNMeanXnParticle_;
        scalarField rhoMMeanXnParticle_;
        scalarField linearKEMeanXnParticle_;
        vectorField momentumMeanXnParticle_;
        List<scalarField> nParcels_;

        // boundary measurements
        List<labelList> boundaryCells_;
        
        List<scalarField> rhoNBF_;
        List<scalarField> rhoMBF_;
        List<scalarField> linearKEBF_;
        List<vectorField> momentumBF_;

        // fields required for decomposition
        volScalarField CLB_;
        volVectorField CLBQ_;
        volTensorField CLBS_;
        volScalarField mu_;
        volScalarField kappa_;
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

    //- Runtime type information
    TypeName("constitutiveLaws");

    //- Construct from components
    constitutiveLaws
    (
        const dictionary& dict,
        const polyMesh& mesh_,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~constitutiveLaws() = default;

    // Member Functions

    void decompose();
 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "uspDynamicLoadBalancingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

