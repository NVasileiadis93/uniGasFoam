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

Class
    Foam::simplifiedChapmanEnskog

Description
    Hybrid decomposition class based on simplified Chapman Enskog criterion

SourceFiles
    simplifiedChapmanEnskog.C

\*---------------------------------------------------------------------------*/

#ifndef simplifiedChapmanEnskog_H
#define simplifiedChapmanEnskog_H

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
                           Class simplifiedChapmanEnskog Declaration
\*---------------------------------------------------------------------------*/

class simplifiedChapmanEnskog
:
    public uspHybridDecomposition
{

    // Private data

        // Properties dictionary
        dictionary propsDict_;

        // Breakdown parameter threshold
        scalar breakdownMax_;
        
        // Time averaging coefficient
        scalar theta_;
        
        // Number of smoothing passes
        label smoothingPasses_;

        label timeSteps_;
        
        scalar timeAvCounter_;

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
        scalarField rhoNMeanXnParticle_;
        scalarField rhoMMeanXnParticle_;
        scalarField linearKEMeanXnParticle_;
        vectorField momentumMeanXnParticle_;
        List<scalarField> nParcels_;
        
        // fields required for decomposition
        volScalarField B_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField translationalT_;
        volVectorField UMean_;
        volVectorField heatFluxVector_;
        volTensorField pressureTensor_;
        volTensorField shearStressTensor_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("simplifiedChapmanEnskog");

    //- Construct from components
    simplifiedChapmanEnskog
    (
        const Time& t,
        const polyMesh& mesh_,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~simplifiedChapmanEnskog() = default;

    // Member Functions

    void decompose();

    void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "uspDynamicLoadBalancingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

