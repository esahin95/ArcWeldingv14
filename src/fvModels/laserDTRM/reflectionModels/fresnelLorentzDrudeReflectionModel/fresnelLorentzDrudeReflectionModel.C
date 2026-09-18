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

#include "fresnelLorentzDrudeReflectionModel.H"
#include "addToRunTimeSelectionTable.H"

#include "fundamentalConstants.H"
#include "universalConstants.H"
#include "electromagneticConstants.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace reflectionModels
    {
        defineTypeNameAndDebug(fresnelLorentzDrude, 0);
        addToRunTimeSelectionTable
        (
            reflectionModel,
            fresnelLorentzDrude,
            dictionary
        );
    }
}

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reflectionModels::fresnelLorentzDrude::fresnelLorentzDrude
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    reflectionModel(dict, mesh),
    N_
    (
        "N",
        dimless/dimVolume,
        dict.lookup<scalar>("N")
    ),
    lambda_
    (
        "lambda",
        dimLength,
        dict.lookup<scalar>("lambda")
    ),
    sigma_
    (
        "sigma",
        dimensionSet(-1,-3,3,0,0,2,0),
        dict.lookup<scalar>("sigma")
    ),
    n_(0),
    k_(0)
{
    // Plasma frequency squared
    const dimensionedScalar omegaPSqr =
        N_*sqr(electromagnetic::e) / electromagnetic::epsilon0 / atomic::me;

    // Damping frequency
    const dimensionedScalar gamma =
        omegaPSqr * electromagnetic::epsilon0 / sigma_;
    Info<<electromagnetic::epsilon0.dimensions()<<endl;

    // Laser angular frequency
    const dimensionedScalar omega =
        mathematical::twoPi * universal::c / lambda_;

    // Relative permittivity
    const scalar epsReal =
        1 - (omegaPSqr / (sqr(gamma)+sqr(omega))).value();
    const scalar epsImag =
        (gamma/omega * omegaPSqr / (sqr(gamma)+sqr(omega))).value();

    // Complex refractive index
    n_ = sqrt(0.5 * (sqrt(sqr(epsReal)+sqr(epsImag)) + epsReal));
    k_ = sqrt(0.5 * (sqrt(sqr(epsReal)+sqr(epsImag)) - epsReal));
    DebugInfo<<"Reflextive index set to n = " << n_
             << " , k = " << k_ << endl;
    DebugInfo<<"Normal incidence reflectivity R = " << rho(1.0) << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::scalar Foam::reflectionModels::fresnelLorentzDrude::rho
(
    const scalar cosTheta
) const
{
    const scalar sinThetaSqr = 1 - sqr(cosTheta);
    const scalar sinThetaTanTheta = sinThetaSqr/cosTheta;

    const scalar nSqr = sqr(n_);
    const scalar kSqr = sqr(k_);

    const scalar x = nSqr - kSqr - sinThetaSqr;
    const scalar aSqr = 0.5 * (sqrt(sqr(x) + 4*nSqr*kSqr) + x);
    const scalar bSqr = 0.5 * (sqrt(sqr(x) + 4*nSqr*kSqr) - x);
    const scalar a = sqrt(aSqr);

    const scalar Rs = (sqr(a-cosTheta) + bSqr) / (sqr(a+cosTheta) + bSqr);
    const scalar Rp = Rs * (sqr(a-sinThetaTanTheta) + bSqr)
                         / (sqr(a+sinThetaTanTheta) + bSqr);

    return 0.5 * (Rs + Rp);
}


// ************************************************************************* //
