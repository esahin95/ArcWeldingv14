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
#include "fvcGrad.H"

#include "physicoChemicalConstants.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace evaporationModels
{
    defineTypeNameAndDebug(Knight, 0);

    addToRunTimeSelectionTable
    (
        evaporationModel,
        Knight,
        dictionary
    );
}
}

using namespace Foam::constant;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evaporationModels::Knight::Knight
(
    const fvMesh& mesh,
    const word& group
)
:
    gasDynamic(mesh, group),

    T0_
    (
        "T0",
        dimTemperature,
        dict_.lookup<scalar>("T0")
    ),

    g0_
    (
        "g0",
        dimless,
        dict_.lookup<scalar>("g0")
    ),

    R0_
    (
        "R0",
        physicoChemical::R/
        dimensionedScalar
        (
            "M",
            dimMass/dimMoles,
            dict_.lookup<scalar>("M0")*1e-3
        )
    ),

    Th_
    (
        "Th",
        dimTemperature,
        dict_.lookup<scalar>("Th")
    ),

    gv_
    (
        "gv",
        dimless,
        dict_.lookup<scalar>("gv")
    )
{
    if (evaporationModel::debug)
    {
        Info<< "R0 = " << R0_ << endl;
    }

    // Update mass transfer rate and recoil pressure
    correct(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::evaporationModels::Knight::correct(const bool relax)
{
    // cache old mass transfer rate
    const volScalarField::Internal mDot0(mDot_);

    // Temporary field data
    tmp<volScalarField::Internal> tcData
    (
        volScalarField::Internal::New("cData", mesh_, dimless)
    );
    volScalarField::Internal& cData = tcData.ref();

    // Modified Mach number
    const volScalarField::Internal m
        = max(min((T_.v() - Tv_)/(Th_ - Tv_), 1.0), 0.0)*sqrt(0.5*gv_);

    cData = 0.5*m*(gv_ - 1.0)/(gv_ + 1.0);
    const volScalarField::Internal sqrtTByTs =
        sqrt(1.0 + mathematical::pi*sqr(cData)) - sqrt(mathematical::pi)*cData;

    cData = sqrt(2.0*Rv_/R0_/g0_)*sqrt(T_.v()/T0_)*sqrtTByTs*m;
    const dimensionedScalar b = (g0_ + 1.0)/4.0;
    const volScalarField::Internal pByP1 =
        1.0 + g0_*cData*(b*cData + sqrt(1.0 + sqr(b*cData)));

    mDot_ = sqrt(2.0/Rv_/T_.v())*m*pByP1/sqrtTByTs*p0_;
    if (relax)
    {
        mDot_ = relax_*mDot_ + (1 - relax_)*mDot0;
    }

    pRec_ = ((1.0 + 2.0*sqr(m))*pByP1 - 1.0)*p0_;

    // Return maximum change
    return gMax(mag(mDot_ - mDot0)().primitiveField());
}


// ************************************************************************* //
