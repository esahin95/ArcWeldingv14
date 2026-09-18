/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "KnightEvaporationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace evaporationModels
{
    defineTypeNameAndDebug(Knight, 0);
    addToRunTimeSelectionTable(evaporationModel, Knight, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evaporationModels::Knight::Knight
(
    const fvMesh& mesh,
    const word& group
)
:
    gasDynamic(mesh, group),
    T0_("T0", dimTemperature, dict_.lookup<scalar>("T0")),
    g0_("g0", dimless, dict_.lookup<scalar>("g0")),
    R0_
    (
        "R0",
        constant::physicoChemical::R
       /dimensionedScalar
        (
            "M",
            dimMass/dimMoles,
            dict_.lookup<scalar>("M0")*1e-3
        )
    ),
    Th_("Th", dimTemperature, dict_.lookup<scalar>("Th")),
    gv_("gv", dimless, dict_.lookup<scalar>("gv"))
{
    if (evaporationModel::debug)
    {
        Info<< "R0 = " << R0_ << endl;
    }

    correct(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::evaporationModels::Knight::correct(const bool relax)
{
    using constant::mathematical::pi;

    const volScalarField::Internal& T = T_.v();
    const volScalarField::Internal mDot0(mDot_);

    // Mach number of the vapour
    const volScalarField::Internal m
    (
        max(min((T - Tv_)/(Th_ - Tv_), 1.0), 0.0)*sqrt(0.5*gv_)
    );

    // Jump conditions across the Knudsen layer
    const volScalarField::Internal c1(0.5*m*(gv_ - 1.0)/(gv_ + 1.0));

    const volScalarField::Internal sqrtTByTs
    (
        sqrt(1.0 + pi*sqr(c1)) - sqrt(pi)*c1
    );

    const volScalarField::Internal c2
    (
        sqrt(2.0*Rv_/R0_/g0_)*sqrt(T/T0_)*sqrtTByTs*m
    );

    const dimensionedScalar b = (g0_ + 1.0)/4.0;

    const volScalarField::Internal pByP1
    (
        1.0 + g0_*c2*(b*c2 + sqrt(1.0 + sqr(b*c2)))
    );

    mDot_ = sqrt(2.0/Rv_/T)*m*pByP1/sqrtTByTs*p0_;

    if (relax)
    {
        mDot_ = relax_*mDot_ + (1 - relax_)*mDot0;
    }

    pRec_ = ((1.0 + 2.0*sqr(m))*pByP1 - 1.0)*p0_;

    return gMax(mag(mDot_.primitiveField() - mDot0.primitiveField()));
}


// ************************************************************************* //
