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
    Foam::uniGasConfiguration

Description
    Abstract class for initialising particles.

SourceFiles
    uniGasConfiguration.C

\*---------------------------------------------------------------------------*/

#ifndef uniGasConfiguration_H
#define uniGasConfiguration_H

#include "Random.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

class uniGasCloud;

/*---------------------------------------------------------------------------*\
                      Class uniGasConfiguration Declaration
\*---------------------------------------------------------------------------*/

class uniGasConfiguration
{
protected:

    // Protected data

        const fvMesh& mesh_;

        uniGasCloud& cloud_;

        const dictionary& uniGasInitialiseDict_;

        Random& rndGen_;

        label nParcelsAdded_;


public:

    //- Runtime type information
    TypeName("uniGasConfiguration");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        uniGasConfiguration,
        dictionary,
        (
            uniGasCloud& cloud,
            const dictionary& dict
        ),
        (cloud, dict)
    );

    //- Construct from components
    uniGasConfiguration(uniGasCloud& cloud, const dictionary& dict);


    //- Selectors
    static autoPtr<uniGasConfiguration> New
    (
        uniGasCloud& cloud,
        const dictionary& dict
    );


    //- Destructor
    virtual ~uniGasConfiguration() = default;


    // Member Functions

        virtual void setInitialConfiguration() = 0;

        virtual label nParcelsAdded() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
