/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
#include "fvcMeshPhi.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::icoMultiPhaseChangeVoF::thermophysicalPredictor()
{
    volScalarField& T = mixture.T();

    const surfaceScalarField rhoPhiCp
    (
        "rhoPhiCp",
        rhoPhi*fvc::interpolate(rhoCp/mixture.rho())
    );


    for(label i=0; i<nThermoCorr_; i++)
    {
        T.storePrevIter();

        fvScalarMatrix TEqn
        (
           fvm::ddt(rhoCp, T) + fvm::div(rhoPhiCp, T)
         - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
         - fvm::laplacian(mixture.kappaEff(momentumTransport.nut()), T)
        ==
          fvModels().source(rhoCp, T)
        );

        // Add phase change contribution
        addSup(TEqn);

        TEqn.relax();

        fvConstraints().constrain(TEqn);

        TEqn.solve();

        fvConstraints().constrain(T);

        // Termination criteria
        if (correctPhaseChange() < thermoTol_ && i > pimple.nCorrNonOrth())
        {
            Info<< "Converged at it = " << i << endl;
            break;
        }
    }

    mixture.correctThermo();
    mixture.correct();
}


// ************************************************************************* //
