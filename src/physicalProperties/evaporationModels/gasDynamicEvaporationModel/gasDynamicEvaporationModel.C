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

#include "gasDynamicEvaporationModel.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcGrad.H"
#include "fvmSup.H"
#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"
#include "fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace evaporationModels
{
    defineTypeNameAndDebug(gasDynamic, 0);

    addToRunTimeSelectionTable
    (
        evaporationModel,
        gasDynamic,
        dictionary
    );
}
}

using namespace Foam::constant;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evaporationModels::gasDynamic::gasDynamic
(
    const fvMesh& mesh,
    const word& group
)
:
    evaporationModel(mesh, group),

    dict_(subDict(dictName_)),

    Lv_
    (
        "Lv",
        dimEnergy/dimMass,
        dict_.lookup<scalar>("L")
    ),

    p0_
    (
        "p0",
        dimPressure,
        dict_.lookup<scalar>("p0")
    ),

    Rv_
    (
        "Rv",
        physicoChemical::R/
        dimensionedScalar
        (
            "M",
            dimMass/dimMoles,
            dict_.lookup<scalar>("Mv")*1e-3
        )
    ),

    Tv_
    (
        "T",
        dimTemperature,
        dict_.lookup<scalar>("Tv")
    ),

    relax_(dict_.lookup<scalar>("relax")),

    pRec_
    (
        IOobject
        (
            IOobject::groupName("pRec", group),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimPressure, 0.0)
    )
{
    if (evaporationModel::debug)
    {
        Info<< "Rv = " << Rv_ << endl;
    }

    // Update mass transfer rate and recoil pressure
    correct(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::evaporationModels::gasDynamic::correct(const bool relax)
{
    // cache old mass transfer rate
    const volScalarField::Internal mDot0(mDot_);

    // Saturated vapor pressure
    const dimensionedScalar LByRTv = Lv_/(Rv_*Tv_);
    const volScalarField::Internal pSat = p0_*exp(LByRTv*(1. - Tv_/T_.v()));

    // Mass transfer rate
    mDot_ = 0.816*pSat/sqrt(mathematical::twoPi*Rv_*T_.v());
    if (relax)
    {
        mDot_ = relax_*mDot_ + (1-relax_)*mDot0;
    }

    // Recoil pressure
    pRec_ = 0.54*pSat;

    // Return maximum change
    return gMax(mag(mDot_ - mDot0)().primitiveField());
}


void Foam::evaporationModels::gasDynamic::addSup(fvMatrix<scalar>& eqn) const
{
    const volScalarField::Internal magGradAlpha = mag(fvc::grad(alpha_)()());

    if (evaporationModel::debug)
    {
        const dimensionedScalar hv = fvc::domainIntegrate
            (
                Lv_*mDot_*magGradAlpha
            );
        Info<< "Total evaporative enthalphy: " << hv << endl;
    }

    // Add latent heat of evaporation
    eqn += Lv_*mDot_*magGradAlpha;
}


void Foam::evaporationModels::gasDynamic::addSup(fvMatrix<vector>& eqn) const
{
    // Add recoil pressure term
    eqn -= pRec_*fvc::grad(alpha_)()();
}

// ************************************************************************* //
