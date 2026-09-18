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

#include "gasDynamicImplicitEvaporationModel.H"
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
    defineTypeNameAndDebug(gasDynamicImplicit, 0);

    addToRunTimeSelectionTable
    (
        evaporationModel,
        gasDynamicImplicit,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evaporationModels::gasDynamicImplicit::gasDynamicImplicit
(
    const fvMesh& mesh,
    const word& group
)
:
    gasDynamic(mesh, group)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::evaporationModels::gasDynamicImplicit::addSup(fvMatrix<scalar>& eqn) const
{
    const volScalarField::Internal& T = T_.v();

    const volScalarField::Internal magGradAlpha = mag(fvc::grad(alpha_)()());
    const fvMatrix<scalar> evapEqn
    (
        fvm::Sp(Lv_*mDot_/T*(Lv_/Rv_/T - 0.5)*magGradAlpha, eqn.psi())
        - Lv_*mDot_*(Lv_/Rv_/T - 1.5)*magGradAlpha
    );

    if (evaporationModel::debug)
    {
        const dimensionedScalar hv = fvc::domainIntegrate(evapEqn&T_);
        Info<< "Total evaporative enthalphy: " << hv << endl;
    }

    eqn += evapEqn;
}


// ************************************************************************* //
