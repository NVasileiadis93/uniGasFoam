/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::uspParcel::constantProperties::constantProperties()
:
    mass_(0),
    d_(0)
{}


inline Foam::uspParcel::constantProperties::constantProperties
(
    const dictionary& dict
)
:
    mass_(dict.get<scalar>("mass")),
    d_(dict.get<scalar>("diameter")),
    rotationalDoF_(dict.get<label>("rotationalDegreesOfFreedom")),
    vibrationalDoF_(dict.get<label>("vibrationalModes")),
    omega_(dict.get<scalar>("omega")),
    alpha_(dict.get<scalar>("alpha")),
    Tref_(dict.get<scalar>("Tref")),
    thetaV_(dict.lookup("characteristicVibrationalTemperature")),
    thetaD_(dict.lookup("dissociationTemperature")),
    Zref_(dict.lookup("Zref")),
    TrefZv_(dict.lookup("referenceTempForZref")),
    charDissQuantumLevel_(dict.lookup("charDissQuantumLevel")),
    ionisationTemperature_(dict.get<scalar>("ionisationTemperature")),
    nElectronicLevels_(dict.get<label>("numberOfElectronicLevels")),
    charge_(dict.get<label>("charge")),
    electronicEnergyList_(dict.lookup("electronicEnergyList")),
    degeneracyList_(dict.lookup("degeneracyList"))
{
    if (thetaV_.size() != vibrationalDoF_)
    {
        FatalErrorInFunction
            << "Number of characteristic vibrational temperatures is "
            << thetaV_.size()
            << ", instead of "
            << vibrationalDoF_ << nl
            << exit(FatalError);
    }

    if (thetaD_.size() != vibrationalDoF_)
    {
        FatalErrorInFunction
            << "Number of characteristic dissociation temperatures is"
            << thetaD_.size()
            << ", instead of "
            << vibrationalDoF_ << nl
            << exit(FatalError);
    }

    if (Zref_.size() != vibrationalDoF_)
    {
        FatalErrorInFunction
            << "Number of reference vibrational relaxation numbers is"
            << Zref_.size()
            << ", instead of "
            << vibrationalDoF_ << nl
            << exit(FatalError);
    }

    if (TrefZv_.size() != vibrationalDoF_)
    {
        FatalErrorInFunction
            << "Number of reference temperature for vibrational relaxation is"
            << TrefZv_.size()
            << ", instead of "
            << vibrationalDoF_ << nl
            << exit(FatalError);
    }

    if (charDissQuantumLevel_.size() != vibrationalDoF_)
    {
        FatalErrorInFunction
            << "Number of characteristic dissociation quantum levels is"
            << charDissQuantumLevel_.size()
            << ", instead of "
            << vibrationalDoF_ << nl
            << exit(FatalError);
    }

    if (degeneracyList_.size() != nElectronicLevels_)
    {
        FatalErrorInFunction
            << "Number of degeneracy levels should be "
            << nElectronicLevels_
            << ", instead of "
            << degeneracyList_.size() << nl
            << exit(FatalError);
    }

    if (electronicEnergyList_.size() != nElectronicLevels_)
    {
        FatalErrorInFunction
            << "Number of electronic energy levels should be "
            << nElectronicLevels_
            << ", instead of "
            << electronicEnergyList_.size() << nl
            << exit(FatalError);
    }

    if (charge_ > 1 || charge_ < -1)
    {
        FatalErrorInFunction
            << "Charge value should be 0 for neutrals, "
            << "1 for ions, or - 1 for electrons"
            << ", instead of "
            << charge_ << nl
            << exit(FatalError);
    }
}


inline Foam::uspParcel::uspParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const vector& U,
    const scalar CWF,
    const scalar RWF,
    const scalar ERot,
    const label ELevel,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label typeId,
    const label newParcel,
    const labelList& vibLevel
)
:
    particle(mesh, coordinates, celli, tetFacei, tetPti),
    U_(U),
    CWF_(CWF),
    RWF_(RWF),
    ERot_(ERot),
    ELevel_(ELevel),
    typeId_(typeId),
    newParcel_(newParcel),
    vibLevel_(vibLevel)
{}


inline Foam::uspParcel::uspParcel
(
    const polyMesh& mesh,
    const vector& position,
    const vector& U,
    const scalar CWF,
    const scalar RWF,
    const scalar ERot,
    const label ELevel,
    const label cellI,
    const label typeId,
    const label newParcel,
    const labelList& vibLevel
)
:
    particle(mesh, position, cellI),
    U_(U),
    CWF_(CWF),
    RWF_(RWF),
    ERot_(ERot),
    ELevel_(ELevel),
    typeId_(typeId),
    newParcel_(newParcel),
    vibLevel_(vibLevel)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

inline Foam::scalar Foam::uspParcel::constantProperties::mass() const
{
    return mass_;
}


inline Foam::scalar Foam::uspParcel::constantProperties::d() const
{
    return d_;
}


inline Foam::scalar Foam::uspParcel::constantProperties::sigmaT() const
{
    return constant::mathematical::pi*d_*d_;
}


inline Foam::label Foam::uspParcel::constantProperties::rotationalDoF() const
{
    return rotationalDoF_;
}


inline Foam::label Foam::uspParcel::constantProperties::vibrationalDoF() const
{
    return vibrationalDoF_;
}


inline Foam::scalar Foam::uspParcel::constantProperties::omega() const
{
    return omega_;
}


inline Foam::scalar Foam::uspParcel::constantProperties::alpha() const
{
    return alpha_;
}


inline Foam::scalar Foam::uspParcel::constantProperties::Tref() const
{
    return Tref_;
}


inline const Foam::scalarList&
Foam::uspParcel::constantProperties::thetaV() const
{
    return thetaV_;
}


inline const Foam::scalarList&
Foam::uspParcel::constantProperties::thetaD() const
{
    return thetaD_;
}


inline const Foam::scalarList&
Foam::uspParcel::constantProperties::Zref() const
{
    return Zref_;
}


inline const Foam::scalarList&
Foam::uspParcel::constantProperties::TrefZv() const
{
    return TrefZv_;
}


inline const Foam::labelList&
Foam::uspParcel::constantProperties::charDissQuantumLevel() const
{
    return charDissQuantumLevel_;
}


inline Foam::scalar
Foam::uspParcel::constantProperties::ionisationTemperature() const
{
    return ionisationTemperature_;
}


inline Foam::label
Foam::uspParcel::constantProperties::nElectronicLevels() const
{
    return nElectronicLevels_;
}


inline Foam::label Foam::uspParcel::constantProperties::charge() const
{
    return charge_;
}


inline const Foam::scalarList&
Foam::uspParcel::constantProperties::electronicEnergyList() const
{
    return electronicEnergyList_;
}


inline const Foam::labelList&
Foam::uspParcel::constantProperties::degeneracyList() const
{
    return degeneracyList_;
}


// * * * * * * * * * * UspParcel Member Functions  * * * * * * * * * * //

inline Foam::label Foam::uspParcel::typeId() const
{
    return typeId_;
}


inline const Foam::vector& Foam::uspParcel::U() const
{
    return U_;
}


inline Foam::scalar Foam::uspParcel::CWF() const
{
    return CWF_;
}


inline Foam::scalar Foam::uspParcel::RWF() const
{
    return RWF_;
}


inline Foam::scalar Foam::uspParcel::ERot() const
{
    return ERot_;
}


inline Foam::label Foam::uspParcel::ELevel() const
{
    return ELevel_;
}


inline Foam::label Foam::uspParcel::newParcel() const
{
    return newParcel_;
}


inline const Foam::labelList& Foam::uspParcel::vibLevel() const
{
    return vibLevel_;
}


inline Foam::label& Foam::uspParcel::typeId()
{
    return typeId_;
}


inline Foam::vector& Foam::uspParcel::U()
{
    return U_;
}


inline Foam::scalar& Foam::uspParcel::CWF()
{
    return CWF_;
}


inline Foam::scalar& Foam::uspParcel::RWF()
{
    return RWF_;
}


inline Foam::scalar& Foam::uspParcel::ERot()
{
    return ERot_;
}


inline Foam::label& Foam::uspParcel::ELevel()
{
    return ELevel_;
}


inline Foam::label& Foam::uspParcel::newParcel()
{
    return newParcel_;
}


inline Foam::labelList& Foam::uspParcel::vibLevel()
{
    return vibLevel_;
}


// ************************************************************************* //
