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
    EllipsoidalStatisticalRelaxation

Description
    Relaxes particles to local Chapman-Enskog distributions

SourceFiles
    EllipsoidalStatisticalRelaxation.C

\*---------------------------------------------------------------------------*/

#ifndef EllipsoidalStatisticalRelaxation_H
#define EllipsoidalStatisticalRelaxation_H

#include "relaxationModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "uspCloud.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        class EllipsoidalStatisticalRelaxation Declaration
\*---------------------------------------------------------------------------*/

class EllipsoidalStatisticalRelaxation
:
    public relaxationModel
{
    // Private data

        label infoCounter_;

        labelList typeIds_;

        scalar rhoNMean_;
        scalar rhoNMeanXnParticle_;
        scalar rhoNMeanInt_;
        scalar molsElec_;
        scalar rhoMMean_;
        scalar rhoMMeanXnParticle_;
        scalar linearKEMean_;
        scalar linearKEMeanXnParticle_;
        scalar rotationalEMean_;
        scalar rotationalDofMean_;
        scalar muu_;
        scalar muv_;
        scalar muw_;
        scalar mvv_;
        scalar mvw_;
        scalar mww_;
        scalar mcc_;
        scalar mccu_;
        scalar mccv_;
        scalar mccw_;
        scalar eu_;
        scalar ev_;
        scalar ew_;
        scalar e_;
        vector momentumMean_;
        vector momentumMeanXnParticle_;
        List<List<scalar>> vibrationalETotal_;
        List<scalar> electronicETotal_;
        List<scalar> nParcels_;
        List<scalar> nParcelsXnParticle_;
        List<scalar> mccSpecies_;
        List<scalar> nGroundElectronicLevel_;
        List<scalar> nFirstElectronicLevel_;

        scalar viscosity_;
        scalar Prandtl_;
        scalar relaxFreq_;
        scalar rhoN_;
        scalar p_;
        scalar translationalT_;
        scalar rotationalT_;
        scalar vibrationalT_;
        scalar electronicT_;
        scalar overallT_;
        vector UMean_;
        tensor pressureTensor_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("EllipsoidalStatisticalRelaxation");

    //- Construct from components
    EllipsoidalStatisticalRelaxation
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~EllipsoidalStatisticalRelaxation() = default;


    // Member Functions

        //- Calculate macroscopic quantities
        virtual void calculateProperties(const label& cell);

        //- Reset macroscopic quantities
        virtual void resetProperties(const label& cell);

        //- Change post-collision velocities to conserve momentum and energy
        virtual void conserveMomentumAndEnergy(const label& cell);

        //- Relax particles
        virtual void relax();

        //- Sample linear velocity from CE distribution
        vector samplePostRelaxationVelocity
        (   
            const scalar& m,
            const scalar& Pr,
            const scalar& p,
            const scalar& T,
            const vector& U,
            const tensor& pT
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

