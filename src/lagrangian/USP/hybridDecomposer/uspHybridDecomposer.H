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

        label timeSteps_;
        label nAvTimeSteps_;

        label decomposeInterval_;
        scalar bMax_;
        scalar Tref_;
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
        List<scalarField> nParcelsXnParticle_;
        List<scalarField> mccSpecies_;
        List<scalarField> nGroundElectronicLevel_;
        List<scalarField> nFirstElectronicLevel_;

        // fields required for decomposition
        volScalarField breakdownParameter_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField translationalT_;
        volScalarField rotationalT_;
        volScalarField vibrationalT_;
        volScalarField electronicT_;
        volVectorField UMean_;
        volVectorField heatFluxVector_;
        volTensorField pressureTensor_;
        volTensorField shearStressTensor_;

    // Private Member Functions

public:

    //- Construct from components
    uspHybridDecomposer
    (
        const dictionary& dict,
        const fvMesh& mesh,
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

