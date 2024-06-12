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
    Foam::cellMeasurements

Description
    This class measures the microscopic properties of particles inside cells. 
    This is required to return internal macroscopic properties.
    
    Reference:
        G. A. Bird. "Visual DSMC Program for Two-Dimensional and Axially
        Symmetric Flows", The DS2V Program User’s Guide, Version 3.8,
        GAB Consulting, Sydney, Australia, October 2006.

SourceFiles
    cellMeasurementsI.H
    cellMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef cellMeasurements_H
#define cellMeasurements_H

#include "fvMesh.H"
#include "uniGasParcel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cellMeasurements Declaration
\*---------------------------------------------------------------------------*/

class cellMeasurements
{
    // Private data

        const fvMesh& mesh_;

        uniGasCloud& cloud_;

        dictionary propsDict_;

        labelList typeIds_;

        List<scalarField> rhoNMean_;

        List<scalarField> rhoNMeanXnParticle_;

        List<scalarField> rhoNMeanInt_;

        List<scalarField> molsElec_;

        List<scalarField> rhoMMean_;

        List<scalarField> rhoMMeanXnParticle_;

        List<scalarField> linearKEMean_;

        List<scalarField> linearKEMeanXnParticle_;

        List<scalarField> rotationalEMean_;

        List<scalarField> rotationalDofMean_;

        List<scalarField> muu_;

        List<scalarField> muv_;

        List<scalarField> muw_;

        List<scalarField> mvv_;

        List<scalarField> mvw_;

        List<scalarField> mww_;

        List<scalarField> mcc_;

        List<scalarField> mccu_;

        List<scalarField> mccv_;

        List<scalarField> mccw_;

        List<scalarField> eRotu_;

        List<scalarField> eRotv_;

        List<scalarField> eRotw_;

        List<scalarField> eRot_;

        List<scalarField> eVibu_;

        List<scalarField> eVibv_;

        List<scalarField> eVibw_;

        List<scalarField> eVib_;

        List<vectorField> momentumMean_;

        List<vectorField> momentumMeanXnParticle_;

        List<List<scalarField>> vibrationalETotal_;

        List<scalarField> electronicETotal_;

        List<scalarField> nParcels_;

        List<scalarField> nParcelsXnParticle_;

        List<scalarField> mccSpecies_;

        List<scalarField> nGroundElectronicLevel_;

        List<scalarField> nFirstElectronicLevel_;

public:

    friend class uniGasParcel;

    // Constructors

        //- Construct from mesh and cloud
        cellMeasurements
        (
            const polyMesh&,
            uniGasCloud& cloud
        );

        //- Construct from mesh, cloud and boolean (uniGasFoam)
        cellMeasurements
        (
            const polyMesh&,
            uniGasCloud& cloud,
            const bool dummy
        );


    //- Destructor
    ~cellMeasurements() = default;

    // Member Functions

        void createFields();

        void calculateFields();

        void clean();

        inline const List<scalarField>& rhoNMean() const;

        inline const List<scalarField>& rhoNMeanXnParticle() const;

        inline const List<scalarField>& rhoNMeanInt() const;

        inline const List<scalarField>& molsElec() const;

        inline const List<scalarField>& rhoMMean() const;

        inline const List<scalarField>& rhoMMeanXnParticle() const;

        inline const List<scalarField>& linearKEMean() const;

        inline const List<scalarField>& linearKEMeanXnParticle() const;

        inline const List<scalarField>& rotationalEMean() const;

        inline const List<scalarField>& rotationalDofMean() const;

        inline const List<scalarField>& muu() const;

        inline const List<scalarField>& muv() const;

        inline const List<scalarField>& muw() const;

        inline const List<scalarField>& mvv() const;

        inline const List<scalarField>& mvw() const;

        inline const List<scalarField>& mww() const;

        inline const List<scalarField>& mcc() const;

        inline const List<scalarField>& mccu() const;

        inline const List<scalarField>& mccv() const;

        inline const List<scalarField>& mccw() const;

        inline const List<scalarField>& eRotu() const;

        inline const List<scalarField>& eRotv() const;

        inline const List<scalarField>& eRotw() const;

        inline const List<scalarField>& eRot() const;

        inline const List<scalarField>& eVibu() const;

        inline const List<scalarField>& eVibv() const;

        inline const List<scalarField>& eVibw() const;

        inline const List<scalarField>& eVib() const;

        inline const List<vectorField>& momentumMean() const;

        inline const List<vectorField>& momentumMeanXnParticle() const;

        inline const List<List<scalarField>>& vibrationalETotal() const;

        inline const List<scalarField>& electronicETotal() const;

        inline const List<scalarField>& nParcels() const;

        inline const List<scalarField>& nParcelsXnParticle() const;

        inline const List<scalarField>& mccSpecies() const;

        inline const List<scalarField>& nGroundElectronicLevel() const;

        inline const List<scalarField>& nFirstElectronicLevel() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
