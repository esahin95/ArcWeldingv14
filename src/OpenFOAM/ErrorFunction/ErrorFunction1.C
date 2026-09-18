/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "ErrorFunction1.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    namespace Function1s
    {
        addScalarFunction1(ErrorFunction);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::ErrorFunction::ErrorFunction
(
    const word& name,
    const unitSets& units,
    const dictionary& dict
)
:
    FieldFunction1<scalar, ErrorFunction>(name),
    Tliq_(dict.lookup<scalar>("Tliq", units.x)),
    Tsol_(dict.lookup<scalar>("Tsol", units.x)),
    Tmid_(0.5 * (Tsol_ + Tliq_)),
    a_(4.0 / (max(small, Tliq_ - Tsol_)))
{
    if (Tliq_ < Tsol_)
    {
        FatalIOErrorInFunction(dict)
            << "Tliq (" << Tliq_ << ") < Tsol (" << Tsol_ << ")"
            << exit(FatalIOError);
    }
}

Foam::Function1s::ErrorFunction::ErrorFunction(const ErrorFunction& errf)
:
    FieldFunction1<scalar, ErrorFunction>(errf),
    Tliq_(errf.Tliq_),
    Tsol_(errf.Tsol_),
    Tmid_(errf.Tmid_),
    a_(errf.a_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1s::ErrorFunction::~ErrorFunction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1s::ErrorFunction::value(const scalar x) const 
{
    scalar y(0.5 * (1.0 - Foam::erf(a_ * (x - Tmid_))));
    return y;
}

Foam::scalar Foam::Function1s::ErrorFunction::integral
(
    const scalar x1, 
    const scalar x2
) const 
{
    NotImplemented;
    return Zero;
}

void Foam::Function1s::ErrorFunction::write
(
    Ostream& os, 
    const unitSets& units
) const 
{
    writeEntry(os, "Tsol", units.x, Tsol_);
    writeEntry(os, "Tliq", units.x, Tliq_);
}