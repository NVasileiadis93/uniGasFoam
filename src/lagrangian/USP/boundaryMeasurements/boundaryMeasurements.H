/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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
    Foam::boundaryMeasurements

Description
    This class measures the microscopic properties of particles that hit boundaries. 
    This is required to return boundary macroscopic properties.
    
    The equations for measuring these properies are given in Bird's 2013 book.
    
    Reference:
        G.A. Bird. "The USP Method", CreateSpace Independent Publishing Platform, 2013.

SourceFiles
    boundaryMeasurementsI.H
    boundaryMeasurements.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryMeasurements_H
#define boundaryMeasurements_H

#include "fvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class uspCloud;
class uspParcel;

/*---------------------------------------------------------------------------*\
                    Class boundaryMeasurements Declaration
\*---------------------------------------------------------------------------*/

class boundaryMeasurements
{
    // Private data

        const fvMesh& mesh_;

        uspCloud& cloud_;

        labelList typeIds_;

        // Outer list is typeId, middle is the patch, inner is the local patch
        // face
        List<List<scalarField>> rhoNIntBF_;
        List<List<scalarField>> rhoNElecBF_;
        List<List<scalarField>> rhoNBF_;
        List<List<scalarField>> rhoMBF_;
        List<List<scalarField>> linearKEBF_;
        List<List<scalarField>> mccSpeciesBF_;
        List<List<vectorField>> momentumBF_;
        List<List<vectorField>> UMeanBF_;
        List<List<scalarField>> rotationalEBF_;
        List<List<scalarField>> rotationalDofBF_;
        List<List<scalarField>> vibrationalEBF_;
        List<List<scalarField>> electronicEBF_;
        List<List<scalarField>> qBF_;
        List<List<vectorField>> fDBF_;


public:

    friend class uspParcel;

    // Constructors

        //- Construct from mesh and cloud
        boundaryMeasurements(const polyMesh&, uspCloud& cloud);

        //- Construct from mesh, cloud and boolean (uspFoam)
        boundaryMeasurements
        (
            const polyMesh&,
            uspCloud& cloud,
            const bool dummy
        );

    //- Destructor
    ~boundaryMeasurements() = default;


    // Member Functions

        void updateFields(uspParcel& p);

        void clean();

        inline const List<List<scalarField>>& rhoNIntBF() const;
        inline const List<List<scalarField>>& rhoNElecBF() const;
        inline const List<List<scalarField>>& rhoNBF() const;
        inline const List<List<scalarField>>& rhoMBF() const;
        inline const List<List<scalarField>>& linearKEBF() const;
        inline const List<List<scalarField>>& mccSpeciesBF() const;
        inline const List<List<vectorField>>& momentumBF() const;
        inline const List<List<vectorField>>& UMeanBF() const;
        inline const List<List<scalarField>>& rotationalEBF() const;
        inline const List<List<scalarField>>& rotationalDofBF() const;
        inline const List<List<scalarField>>& vibrationalEBF() const;
        inline const List<List<scalarField>>& electronicEBF() const;
        inline const List<List<scalarField>>& qBF() const;
        inline const List<List<vectorField>>& fDBF() const;

        inline List<List<scalarField>>& rhoNIntBF();
        inline List<List<scalarField>>& rhoNElecBF();
        inline List<List<scalarField>>& rhoNBF();
        inline List<List<scalarField>>& rhoMBF();
        inline List<List<scalarField>>& linearKEBF();
        inline List<List<scalarField>>& mccSpeciesBF();
        inline List<List<vectorField>>& momentumBF();
        inline List<List<vectorField>>& UMeanBF();
        inline List<List<scalarField>>& rotationalEBF();
        inline List<List<scalarField>>& rotationalDofBF();
        inline List<List<scalarField>>& vibrationalEBF();
        inline List<List<scalarField>>& electronicEBF();
        inline List<List<scalarField>>& qBF();
        inline List<List<vectorField>>& fDBF();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boundaryMeasurementsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
