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
Description
\*----------------------------------------------------------------------------*/

#include "uniGasDynamicLoadBalancing.H"
#include "uniGasCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
uniGasDynamicLoadBalancing::uniGasDynamicLoadBalancing
(
    const Time& t,
    uniGasCloud& cloud
)
:
    time_(t),
    cloud_(cloud)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uniGasDynamicLoadBalancing::calculate
(
)
{

    if ( Pstream::parRun() )
    {
            
        scalar nGlobalParticles = cloud_.size();
        Foam::reduce(nGlobalParticles, sumOp<scalar>());
        
        scalar idealNParticles = scalar(nGlobalParticles)/scalar(Pstream::nProcs());
        
        scalar nParticles = cloud_.size();
        scalar localImbalance = mag(nParticles - idealNParticles);
        Foam::reduce(localImbalance, maxOp<scalar>());
        scalar maxImbalance = localImbalance/idealNParticles;
        Info << "    Maximum imbalance               = " << 100*maxImbalance << endl;
    }

}

}  // End namespace Foam

// ************************************************************************* //