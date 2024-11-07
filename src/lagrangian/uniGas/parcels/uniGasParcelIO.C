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

\*---------------------------------------------------------------------------*/

#include "uniGasParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "uniGasCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniGasParcel::uniGasParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat),
    U_(Zero),
    CWF_(1.0),
    RWF_(1.0),
    ERot_(0.0),
    ELevel_(0),
    typeId_(-1),
    newParcel_(0),
    vibLevel_(0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> U_;
            CWF_ = readScalar(is);
            RWF_ = readScalar(is);
            ERot_ = readScalar(is);
            ELevel_ = readLabel(is);
            typeId_ = readLabel(is);
            newParcel_ = readLabel(is);
            is >> vibLevel_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&U_),
                sizeof(U_)
              + sizeof(CWF_)
              + sizeof(RWF_)
              + sizeof(ERot_)
              + sizeof(ELevel_)
              + sizeof(typeId_)
              + sizeof(newParcel_)
            );
            is >> vibLevel_;
        }
    }

    // Check state of Istream
    is.check(FUNCTION_NAME);
}


void Foam::uniGasParcel::readFields(Cloud<uniGasParcel>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<vector>U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<scalar>CWF(c.fieldIOobject("cellWeight", IOobject::MUST_READ));
    c.checkFieldIOobject(c, CWF);

    IOField<scalar>RWF(c.fieldIOobject("radialWeight", IOobject::MUST_READ));
    c.checkFieldIOobject(c, RWF);

    IOField<scalar>ERot(c.fieldIOobject("ERot", IOobject::MUST_READ));
    c.checkFieldIOobject(c, ERot);

    IOField<label> ELevel(c.fieldIOobject("ELevel", IOobject::MUST_READ));
    c.checkFieldIOobject(c, ELevel);

    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeId);

    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::MUST_READ));
    c.checkFieldIOobject(c, newParcel);

    IOField<labelField> vibLevel(c.fieldIOobject("vibLevel",
                                                 IOobject::MUST_READ));
    c.checkFieldIOobject(c, vibLevel);

    label i = 0;
    forAllIter(uniGasCloud, c, iter)
    {
        uniGasParcel& p = iter();

        p.U_ = U[i];
        p.CWF_ = CWF[i];
        p.RWF_ = RWF[i];
        p.ERot_ = ERot[i];
        p.ELevel_ = ELevel[i];
        p.typeId_ = typeId[i];
        p.newParcel_ = newParcel[i];
        p.vibLevel_ = vibLevel[i];
        i++;
    }
}


void Foam::uniGasParcel::writeFields(const Cloud<uniGasParcel>& c)
{
    particle::writeFields(c);

    label np = c.size();

    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> CWF(c.fieldIOobject("cellWeight", IOobject::NO_READ), np);
    IOField<scalar> RWF(c.fieldIOobject("radialWeight", IOobject::NO_READ), np);
    IOField<scalar> ERot(c.fieldIOobject("ERot", IOobject::NO_READ), np);
    IOField<label> ELevel(c.fieldIOobject("ELevel", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::NO_READ), np);
    IOField<labelField> vibLevel(c.fieldIOobject("vibLevel", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(uniGasCloud, c, iter)
    {
        const uniGasParcel& p = iter();

        U[i] = p.U();
        CWF[i] = p.CWF();
        RWF[i] = p.RWF();
        ERot[i] = p.ERot();
        ELevel[i] = p.ELevel();
        typeId[i] = p.typeId();
        newParcel[i] = p.newParcel();
        vibLevel[i] = p.vibLevel();
        i++;
    }

    U.write();
    CWF.write();
    RWF.write();
    ERot.write();
    vibLevel.write();
    ELevel.write();
    typeId.write();
    newParcel.write();
    vibLevel.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const uniGasParcel& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.U()
            << token::SPACE << p.CWF()
            << token::SPACE << p.RWF()
            << token::SPACE << p.ERot()
            << token::SPACE << p.ELevel()
            << token::SPACE << p.typeId()
            << token::SPACE << p.newParcel()
            << token::SPACE << p.vibLevel();
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.U_),
            sizeof(p.U())
          + sizeof(p.CWF())
          + sizeof(p.RWF())
          + sizeof(p.ERot())
          + sizeof(p.ELevel())
          + sizeof(p.typeId())
          + sizeof(p.newParcel())
        );
        os << p.vibLevel();
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
