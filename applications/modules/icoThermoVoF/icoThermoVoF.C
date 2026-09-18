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

#include "icoThermoVoF.H"
#include "localEulerDdtScheme.H"
#include "fvCorrectPhi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(icoThermoVoF, 0);
    addToRunTimeSelectionTable(solver, icoThermoVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::icoThermoVoF::icoThermoVoF(fvMesh& mesh)
:
    twoPhaseVoFSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new compressibleTwoPhaseVoFMixture(mesh))
    ),

    mixture_
    (
        refCast<compressibleTwoPhaseVoFMixture>(twoPhaseVoFSolver::mixture)
    ),

    p(mixture_.p()),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict(),
        false
    ),

    alphaRhoPhi1
    (
        IOobject::groupName("alphaRhoPhi", alpha1.group()),
        fvc::interpolate(mixture_.thermo1().rho())*alphaPhi1
    ),

    alphaRhoPhi2
    (
        IOobject::groupName("alphaRhoPhi", alpha2.group()),
        fvc::interpolate(mixture_.thermo2().rho())*alphaPhi2
    ),

    rhoCp
    (
        "rhoCp",
        (
            mixture_.thermo1().rho()*alpha1*mixture_.thermo1().Cp()
          + mixture_.thermo2().rho()*alpha2*mixture_.thermo2().Cp()
        )
    ),

    Cp
    (
        IOobject
        (
            "Cp",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rhoCp/rho
    ),

    rhoPhiCp
    (
        "rhoPhiCp",
        rhoPhi * fvc::interpolate(Cp)
    ),

    momentumTransport
    (
        rho,
        U,
        phi,
        rhoPhi,
        alphaPhi1,
        alphaPhi2,
        alphaRhoPhi1,
        alphaRhoPhi2,
        mixture_
    ),

    thermophysicalTransport(momentumTransport)
{
    if (!mixture_.incompressible())
    {
        FatalErrorInFunction
                << "At least one phase is compressible!"
                << exit(FatalError);
    }

    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }

    if (!runTime.restart() || !divergent())
    {
        correctUphiBCs(U_, phi_, true);

        fv::correctPhi
        (
            phi_,
            U,
            p_rgh,
            rAU,
            autoPtr<volScalarField>(),
            pressureReference(),
            pimple
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::icoThermoVoF::~icoThermoVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::icoThermoVoF::maxDeltaT() const
{

    // Diffusion Courant number
    const scalarField sumDif
    (
        fvc::surfaceSum
        (
            mesh.magSf()
            * fvc::interpolate(thermophysicalTransport.kappaEff())
            * mesh.surfaceInterpolation::deltaCoeffs()
        )().primitiveField()
    );

    const scalarField sumDiv
    (
        (
            mesh.V() * rho() * Cp()
        )().primitiveField()
    );

    const scalar diCoNum =
        gMax(sumDif / sumDiv) * runTime.deltaTValue();

    const scalar meanDiCoNum =
        gSum(sumDif) / gSum(sumDiv) * runTime.deltaTValue();

    Info<< "Diffusion Courant Number mean: " << meanDiCoNum
        << " max: " << diCoNum << endl;

    // Recompute maximum time step
    const scalar maxDiCo =
        runTime.controlDict().lookup<scalar>("maxDiCo");

    scalar deltaT = twoPhaseVoFSolver::maxDeltaT();

    if (diCoNum > small)
    {
        deltaT = min(deltaT, maxDiCo/diCoNum*runTime.deltaTValue());
    }

    return deltaT;
}

void Foam::solvers::icoThermoVoF::prePredictor()
{
    twoPhaseVoFSolver::prePredictor();

    const volScalarField& rho1 = mixture_.thermo1().rho();
    const volScalarField& rho2 = mixture_.thermo2().rho();

    const volScalarField& Cp1 = mixture_.thermo1().Cp();
    const volScalarField& Cp2 = mixture_.thermo2().Cp();

    // Phase mass fluxes
    alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;
    alphaRhoPhi2 = fvc::interpolate(rho2)*alphaPhi2;

    // Mass flux
    rhoPhi = alphaRhoPhi1 + alphaRhoPhi2;

    // Heat capacity
    rhoCp = alpha1*rho1*Cp1 + alpha2*rho2*Cp2;
    Cp = rhoCp/rho;

    // Heat capacity flux
    rhoPhiCp = alphaRhoPhi1*fvc::interpolate(Cp1)
             + alphaRhoPhi2*fvc::interpolate(Cp2);
}


void Foam::solvers::icoThermoVoF::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::icoThermoVoF::thermophysicalTransportPredictor()
{
    thermophysicalTransport.predict();
}


void Foam::solvers::icoThermoVoF::pressureCorrector()
{
    incompressiblePressureCorrector(p);
}


void Foam::solvers::icoThermoVoF::momentumTransportCorrector()
{
    momentumTransport.correct();
}


void Foam::solvers::icoThermoVoF::thermophysicalTransportCorrector()
{
    thermophysicalTransport.correct();
}


// ************************************************************************* //
