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
    Foam::uspFaceTracker

Description

SourceFiles
    uspFaceTrackerI.H
    uspFaceTracker.C

\*---------------------------------------------------------------------------*/

#ifndef uspFaceTrackerFO_H
#define uspFaceTrackerFO_H

#include "vector.H"
#include "fvMesh.H"
#include "uspParcel.H"
#include "surfaceFields.H"
#include "CloudFunctionObject.H"
#include "CloudWithModels.H"
#include "uspCloud.H"

#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class uspFaceTracker Declaration
\*---------------------------------------------------------------------------*/

class uspFaceTrackerFO
:
    public CloudFunctionObject<uspCloud>
{
    // Private data

        List<scalarField> parcelIdFlux_;
        List<scalarField> massIdFlux_;


public:

    TypeName("uspFaceTrackerFO");

    // Constructors

        //- Construct from dictionary
        uspFaceTrackerFO
        (
            const dictionary& dict,
            uspCloud& owner,
            const word& modelName
        );

        //- Construct copy
        uspFaceTrackerFO(const uspFaceTrackerFO& ft);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<uspCloud>> clone() const
        {
            return autoPtr<CloudFunctionObject<uspCloud>>
            (
                new uspFaceTrackerFO(*this)
            );
        }

    //- Destructor
    virtual ~uspFaceTrackerFO() = default;


    // Member Functions

        //- Post-evolve hook
        virtual void postEvolve();

        //- Post-face hook
        virtual void postFace(const uspParcel& p, bool& keepParticle);

        inline const List<scalarField>& parcelIdFlux() const;

        inline const List<scalarField>& massIdFlux() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uspFaceTrackerFOI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
