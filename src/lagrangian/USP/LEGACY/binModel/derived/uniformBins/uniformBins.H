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


Description

SourceFiles
    uniformBinsI.H
    uniformBins.C
    uniformBinsIO.C

\*---------------------------------------------------------------------------*/

#ifndef uniformBins_H
#define uniformBins_H


#include "binModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class uniformBins Declaration
\*---------------------------------------------------------------------------*/

class uniformBins
:
    public binModel
{
private:

    // Private data

        dictionary propsDict_;

        vector startPoint_;
        vector endPoint_;
        vector unitVector_;
        scalar rSEMag_;


        label nBins_;
        scalar binWidth_;
        scalar area_;


public:

    //- Runtime type information
    TypeName("uniformBins");

    //- Construct from components
    uniformBins
    (
        const polyMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    ~uniformBins();


    // Member Functions

        label isPointWithinBin
        (
            const vector& rI,
            const label cellI
        );

        scalarField binPositions();

        vectorField bins();

        label nBins() const;

        scalar binVolume(const label n);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
