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
    uspVelocityDistributionZone

Description


SourceFiles
    uspVelocityDistributionZoneI.H
    uspVelocityDistributionZone.C
    uspVelocityDistributionZoneIO.C

\*---------------------------------------------------------------------------*/

#ifndef uspVelocityDistributionZone_H
#define uspVelocityDistributionZone_H

#include "uspField.H"
#include "vector.H"
#include "uspCloud.H"
#include "distribution.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class uspVelocityDistributionZone Declaration
\*---------------------------------------------------------------------------*/

class uspVelocityDistributionZone
:
	public uspField
{
    // Private data

        dictionary propsDict_;

        word fieldName_;

        labelField typeIds_;

        dictionary timeDictVel_;
        timeDataMeas timeVel_;

        word regionName_;
        label regionId_;

        vector UMean_;
        vector Ucollected_;
        label nParcels_;

        scalar binWidth_;
        distribution distrX_;
        distribution distrY_;
        distribution distrZ_;


public:

    //- Runtime type information
    TypeName("uspVelocityDistributionZone");

    //- Construct from components
    uspVelocityDistributionZone
    (
        const Time& t,
        const polyMesh& mesh,
        uspCloud& cloud,
        const dictionary& dict
    );

    //- Destructor
    virtual ~uspVelocityDistributionZone() = default;


    // Member Functions

        //- Create initial field
        virtual void createField();

        //- Calculate field
        virtual void calculateField();

        //- Write field
        virtual void writeField();

        virtual void updateProperties(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
