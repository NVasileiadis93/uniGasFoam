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
    moleculeElectronDissociationIonisation

Description


SourceFiles
    moleculeElectronDissociationIonisation.C


\*---------------------------------------------------------------------------*/

#ifndef moleculeElectronDissociationIonisation_H
#define moleculeElectronDissociationIonisation_H

#include "uniGasReaction.H"
#include "vector.H"
#include "uniGasCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
           Class moleculeElectronDissociationIonisation Declaration
\*---------------------------------------------------------------------------*/

class moleculeElectronDissociationIonisation
:
    public uniGasReaction
{
    // Private data

        labelPair productIdsDiss_;
        labelPair productIdsIon_;

        scalar heatOfReactionDiss_;
        scalar heatOfReactionIon_;

        label nDissociationReactions_;
        label nDissociationReactionsPerTimeStep_;
        label nIonisationReactions_;
        label nIonisationReactionsPerTimeStep_;


    // Private Member Functions

        void setProperties();
        
        void react
        (
            uniGasParcel& p,
            uniGasParcel& q
        );


public:

    //- Runtime type information
    TypeName("moleculeElectronDissociationIonisation");

    // Constructors

        //- Construct from components
        moleculeElectronDissociationIonisation
        (
            const Time& t,
            uniGasCloud& cloud,
            const dictionary& dict
        );

        //- Destructor
        virtual ~moleculeElectronDissociationIonisation() = default;


    // Member Functions

        //- Initial configuration - read in properties from chemReactDict
        void initialConfiguration();

        //- Apply a reaction between parcel p and q
        void reaction(uniGasParcel& p, uniGasParcel& q);

        //- Return true if molecules belong to this reaction model
        bool tryReactMolecules(label typeIdP, const label typeIdQ) const;

        bool outputResults(const label counterIndex);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
