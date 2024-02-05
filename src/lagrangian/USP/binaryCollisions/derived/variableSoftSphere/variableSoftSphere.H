/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.s

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::variableSoftSphere

Description
    This class provides functions to use the variable soft sphere collision
    model. In this model, the effective diameter that allows the collision
    cross-section to be caculated is calculated from

        \f[
            d_{pq} = \left(d_{ref}\right)_{pq} \left[ \left{ 2 k \left(T_{ref}\right)_{pq}/\left( m_r c_r^2\right)\right}^{\omega - \frac{1}{2}}/Gamma\left(5/2 - \omega_{pq} \right) \right]^{\frac{1}{2}}
        \f]

    where
    \vartable
        d_{pq}                    | Effective diameter of collision pair [m]
        \left(d_{ref}\right)_{pq} | Diameter of collision pair at reference temperature [m]
        k                         | Boltzmann constant [m^2 kg s^{-2} K^{-1}]
        T_{ref}                   | Reference temperature [K]
        m_{r}                     | Reduced mass of collision pair [kg]
        c_{r}                     | Relative speed of collision pair [m/s]
        \omega                    | Viscosity index
    \endvartable

    A random elevation and azimuth angle are then generated to calculate the
    post-collision velocities

        \f[
            \cos{\chi} = 2 R_f^{\frac{1}{\alpha}} - 1,
        \f]
        \f[
            \theta = 2 \pi R_f
        \f]

    where R_f is a uniformly distributed random number between 0 and 1.

    Reference:
        G.A. Bird,  "Molecular Gas Dynamics and the Direct Simulation of Gas
        Flows", Oxford Science Publications, 1994 (pages 41, 97, 222, 426, 428).

Usage
    \table
        Property                 | Description                   | Required  | Default
        Tref                     | Reference temperature         | yes       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    binaryCollisionModel            variableSoftSphere;

    collisionProperties
    {
        Tref                        273;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef variableSoftSphere_H
#define variableSoftSphere_H

#include "binaryCollisionModel.H"
#include "uspCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class variableSoftSphere Declaration
\*---------------------------------------------------------------------------*/

class variableSoftSphere
:
    public binaryCollisionModel
{
    // Private data

        //- The properties dictionary
        const dictionary propertiesDict_;

        //- Reference temperature
        const scalar Tref_;


public:

    //- Runtime type information
    TypeName("variableSoftSphere");

    //- Construct from dictionary
    variableSoftSphere
    (
        const dictionary& dict,
        uspCloud& cloud
    );

    //- Destructor
    virtual ~variableSoftSphere() = default;


    // Member Functions

        //- Flag to indicate whether model activates collision model
        virtual bool active() const;

        //- Return the collision cross section * relative velocity product
        virtual scalar sigmaTcR
        (
            const uspParcel& pP,
            const uspParcel& pQ
        ) const;

        //- Apply collision
        virtual void collide
        (
            uspParcel& pP,
            uspParcel& pQ,
            label& cellI
        );

        virtual const dictionary& propertiesDict() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "variableSoftSphere.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
