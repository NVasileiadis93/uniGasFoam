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
    unifiedStochasticParticleSBGK

Description
    Collides particles to local unified SBGK distribution
    
    Use in uniGasProperties file.

Usage
    Example of bgk collision specification

    \verbatim
    bgkCollisionModel    unifiedStochasticParticleSBGK;
    collisionProperties
    {
        Tref                273;
        macroInterpolation  true;
        theta               0.5;
    }
    /endverbatim
    
SourceFiles
    unifiedStochasticParticleSBGK.C

\*---------------------------------------------------------------------------*/

#ifndef unifiedStochasticParticleSBGK_H
#define unifiedStochasticParticleSBGK_H

#include "bgkCollisionModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "uniGasCloud.H"
#include "interpolationCellPoint.H"
#include "interpolationCellPointFace.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        class unifiedStochasticParticleSBGK Declaration
\*---------------------------------------------------------------------------*/

class unifiedStochasticParticleSBGK
:
    public bgkCollisionModel
{
    // Private data
    
        const dictionary propertiesDict_;
        const scalar Tref_;
        const scalar theta_;
        const bool macroInterpolation_;
        label infoCounter_;
        labelList typeIds_;
        scalar shufflePasses_;
        scalar maxProbResetValue_;

        boolField performCollision_;
        boolField maxProbReset_;
        scalarField maxProb_;
        scalarField rhoNMean_;
        scalarField rhoMMean_;
        scalarField linearKEMean_;
        vectorField momentumMean_;
        scalarField rotationalEMean_;
        scalarField rotationalDofMean_;
        scalarField rhoNMeanXnParticle_;
        scalarField rhoMMeanXnParticle_;
        vectorField momentumMeanXnParticle_ ;
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
        scalarField eRotu_;
        scalarField eRotv_;
        scalarField eRotw_;
        scalarField eRot_;
        scalarField eVibu_;
        scalarField eVibv_;
        scalarField eVibw_;
        scalarField eVib_;
        scalarField rhoNMeanInt_;
        scalarField molsElec_;
        List<scalarField> nParcels_;
        List<scalarField> nParcelsXnParticle_;
        List<scalarField> mccSpecies_;
        List<scalarField> nGroundElectronicLevel_;
        List<scalarField> nFirstElectronicLevel_;
        List<scalarField> electronicETotal_;
        List<List<scalarField>> vibrationalETotal_;

        volScalarField Prandtl_;
        volScalarField collFreq_;
        volScalarField rhoN_;
        volScalarField p_;
        volScalarField translationalT_;
        volScalarField rotationalT_;
        volScalarField vibrationalT_;
        volScalarField electronicT_;
        volScalarField overallT_;
        volVectorField UMean_;
        volVectorField heatFluxVector_;
        volVectorField heatFluxVectorPrevious_;
        volTensorField pressureTensor_;
        volTensorField shearStressTensor_;
        volTensorField shearStressTensorPrevious_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("unifiedStochasticParticleSBGK");

    //- Construct from componentsonst label& cell
    unifiedStochasticParticleSBGK
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uniGasCloud& cloud
    );

    //- Destructor
    virtual ~unifiedStochasticParticleSBGK() = default;

    // Member Functions

        //- Return the coefficients dictionary
        virtual const dictionary& propertiesDict() const;

        //- Calculate macroscopic quantities
        virtual void calculateProperties();

        //- Reset macroscopic quantities
        virtual void resetProperties();

        //- Change post-collision velocities to conserve momentum and energy
        virtual void conserveMomentumAndEnergy(const label& cell);

        //- Collide particles
        virtual void collide();

        //- Sample linear velocity from CE distribution
        vector samplePostCollisionVelocity
        (   
            const label& cell,
            const scalar& m,
            const scalar& Pr,
            const scalar& rF,
            const scalar& p,
            const scalar& T,
            const vector& U,
            const vector& q,
            const tensor& s
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

