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
    USBGK

Description
    Relaxes particles to local unified Shakhov distribution
    
    Use in uspProperties file.

Usage
    Example of relaxation specification

    \verbatim
    relaxationCollisionModel    USBGK;
    collisionProperties
    {
        Tref                273;
        macroInterpolation  true;
        theta               0.5;
    }
    /endverbatim
    
SourceFiles
    USBGK.C

\*---------------------------------------------------------------------------*/

#ifndef USBGK_H
#define USBGK_H

#include "relaxationCollisionModel.H"
#include "vector.H"
#include "polyMesh.H"
#include "uspCloud.H"
#include "interpolationCellPoint.H"
#include "interpolationCellPointFace.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        class USBGK Declaration
\*---------------------------------------------------------------------------*/

class USBGK
:
    public relaxationCollisionModel
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

        boolField performRelaxation_;
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
        scalarField eu_;
        scalarField ev_;
        scalarField ew_;
        scalarField e_;
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
        volScalarField relaxFreq_;
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
    TypeName("USBGK");

    //- Construct from componentsonst label& cell
    USBGK
    (
        const dictionary& dict,
        const polyMesh& mesh,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~USBGK() = default;

    // Member Functions

        //- Return the coefficients dictionary
        virtual const dictionary& propertiesDict() const;

        //- Calculate macroscopic quantities
        virtual void calculateProperties();

        //- Reset macroscopic quantities
        virtual void resetProperties();

        //- Change post-collision velocities to conserve momentum and energy
        virtual void conserveMomentumAndEnergy(const label& cell);

        //- Relax particles
        virtual void relax();

        //- Sample linear velocity from CE distribution
        vector samplePostRelaxationVelocity
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

