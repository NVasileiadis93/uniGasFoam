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
    Foam::uspConfiguration

Description
    Abstract class for initialising particles.

SourceFiles
    uspConfiguration.C

\*---------------------------------------------------------------------------*/

#ifndef uspConfiguration_H
#define uspConfiguration_H

#include "Random.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

class uspCloud;

/*---------------------------------------------------------------------------*\
                      Class uspConfiguration Declaration
\*---------------------------------------------------------------------------*/

class uspConfiguration
{
protected:

    // Protected data

        const fvMesh& mesh_;

        uspCloud& cloud_;

        const dictionary& uspInitialiseDict_;

        Random& rndGen_;

        label nParcelsAdded_;


public:

    //- Runtime type information
    TypeName("uspConfiguration");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        uspConfiguration,
        dictionary,
        (
            uspCloud& cloud,
            const dictionary& dict
        ),
        (cloud, dict)
    );

    //- Construct from components
    uspConfiguration(uspCloud& cloud, const dictionary& dict);


    //- Selectors
    static autoPtr<uspConfiguration> New
    (
        uspCloud& cloud,
        const dictionary& dict
    );


    //- Destructor
    virtual ~uspConfiguration() = default;


    // Member Functions

        virtual void setInitialConfiguration() = 0;

        virtual label nParcelsAdded() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
