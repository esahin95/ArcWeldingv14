/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "icoMultiPhaseChangeVoF.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(icoMultiPhaseChangeVoF, 0);
    addToRunTimeSelectionTable(solver, icoMultiPhaseChangeVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::icoMultiPhaseChangeVoF::icoMultiPhaseChangeVoF(fvMesh& mesh)
:
    icoMulticomponentVoF(mesh),

    solidificationModels_(phases.size()),

    evaporationModels_(phases.size()),

    nThermoCorr_
    (
        pimple.dict().lookupOrDefault<label>("nThermoCorr", 20)
    ),

    thermoTol_
    (
        pimple.dict().lookupOrDefault<scalar>("thermoTol", 1e-6)
    )
{
    forAll(phases, phasei)
    {
        solidificationModels_.set
        (
            phasei,
            solidificationModel::New(mesh, phases[phasei].name())
        );

        evaporationModels_.set
        (
            phasei,
            evaporationModel::New(mesh, phases[phasei].name())
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::icoMultiPhaseChangeVoF::~icoMultiPhaseChangeVoF()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::solvers::icoMultiPhaseChangeVoF::correctPhaseChange()
{
    const volScalarField& T = mixture.T();

    // Change of the temperature and of the phase change models
    const scalar resT =
        gMax(mag(T.prevIter().primitiveField() - T.primitiveField()));

    scalar resS = 0;
    forAll(solidificationModels_, phasei)
    {
        resS = max(resS, solidificationModels_[phasei].correct());
    }

    scalar resE = 0;
    forAll(evaporationModels_, phasei)
    {
        resE = max(resE, evaporationModels_[phasei].correct());
    }

    Info<< "resT = " << resT << " , "
        << "resS = " << resS << " , "
        << "resE = " << resE << endl;

    return max(resS, resE);
}


// ************************************************************************* //
