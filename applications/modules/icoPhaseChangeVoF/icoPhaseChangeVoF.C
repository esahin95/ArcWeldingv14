/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "icoPhaseChangeVoF.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(icoPhaseChangeVoF, 0);
    addToRunTimeSelectionTable(solver, icoPhaseChangeVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::icoPhaseChangeVoF::icoPhaseChangeVoF(fvMesh& mesh)
:
    icoThermoVoF(mesh),

    solidificationModel_
    (
        solidificationModel::New(mesh, mixture_.phase1Name())
    ),

    evaporationModel_
    (
        evaporationModel::New(mesh, mixture_.phase1Name())
    ),

    nThermoCorr_
    (
        pimple.dict().lookupOrDefault<label>("nThermoCorr", 20)
    ),

    thermoTol_
    (
        pimple.dict().lookupOrDefault<scalar>("thermoTol", 1e-6)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::icoPhaseChangeVoF::~icoPhaseChangeVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::icoPhaseChangeVoF::correctPhaseChange()
{
    const volScalarField& T = mixture_.T();

    // Change of the temperature and of the phase change models
    const scalar resT =
        gMax(mag(T.prevIter().primitiveField() - T.primitiveField()));
    const scalar resS = solidificationModel_->correct();
    const scalar resE = evaporationModel_->correct();

    Info<< "resT = " << resT << " , "
        << "resS = " << resS << " , "
        << "resE = " << resE << endl;

    return max(resS, resE);
}


// ************************************************************************* //
