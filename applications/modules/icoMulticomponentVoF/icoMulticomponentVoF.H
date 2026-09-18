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

Class
    Foam::solvers::icoMulticomponentVoF

Description
    Solver module for the solution of multiple incompressible, non-isothermal
    immiscible fluids using a VOF (volume of fluid) phase-fraction based
    interface capturing approach, with optional mesh motion and mesh topology
    changes including adaptive re-meshing. Selection of phases can be set to
    be missible with zero surface tension and interface compression, but
    non-zero diffusivity.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    A mixture approach for momentum transport is provided in which a single
    laminar, RAS or LES model is selected to model the momentum stress.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Energy equation is based on icoReactingMultiphaseInterFoam with mixture
    averaged binary mass diffusion coefficients and diffusion fluxes from
    generalized Fick's law. Alpha predictor and pressure corrector are taken
    from incompressibleMultiphaseVoF.

SourceFiles
    icoMulticomponentVoF.C

See also
    Foam::solvers::multiphaseVoFSolver

\*---------------------------------------------------------------------------*/

#ifndef icoMulticomponentVoF_H
#define icoMulticomponentVoF_H

#include "multiphaseVoFSolver.H"
#include "multicomponentVoFMixture.H"
#include "compressibleMomentumTransportModel.H"
#include "pressureReference.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class icoMulticomponentVoF Declaration
\*---------------------------------------------------------------------------*/

class icoMulticomponentVoF
:
    public multiphaseVoFSolver
{

protected:

    // Phase properties

        //- The compressible two-phase mixture
        multicomponentVoFMixture& mixture;

        //- Reference to the phases
        UPtrListDictionary<compressibleVoFphase>& phases;


    // Thermophysical properties

        //- Reference to the mixture static pressure field
        volScalarField& p;

        //- Heat capacities
        const volScalarField& rhoCp;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference_;


    // Kinematic properties

        //- Phase mass fluxes
        PtrList<surfaceScalarField> alphaRhoPhis_;


    // Momentum transport

        //- Momentum transport model pointer
        autoPtr<compressible::momentumTransportModel> momentumTransport_;

        //- Momentum transport model
        compressible::momentumTransportModel& momentumTransport;


    // Cached temporary fields

        //- Continuity error
        tmp<volScalarField> contErr;


    // Protected Member Functions

        //- Return the pressure reference
        virtual const Foam::pressureReference& pressureReference() const
        {
            return pressureReference_;
        }

        //- The flow is incompressible if all phases are incompressible
        virtual bool incompressible() const
        {
            return mixture.incompressible();
        }

        //- The flow is divergent if it is not incompressible
        //  Mass sources are not currently supported
        virtual bool divergent() const
        {
            return !incompressible();
        }

        //- Return the mixture compressibility/density
        //  Not required for incompressible fluids
        virtual tmp<volScalarField> psiByRho() const
        {
            return tmp<volScalarField>(nullptr);
        }

        //- Return the momentum equation stress term
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U);

        //- Solve for the phase-fractions
        virtual void alphaPredictor();

        //- Return the interface surface tension force for the momentum equation
        virtual tmp<surfaceScalarField> surfaceTensionForce() const;


private:

    // Private Member Functions

        //- Solve for the phase-fractions
        void alphaSolve();


public:

    //- Runtime type information
    TypeName("icoMulticomponentVoF");


    // Constructors

        //- Construct from region mesh
        icoMulticomponentVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        icoMulticomponentVoF
        (
            const icoMulticomponentVoF&
        ) = delete;


    //- Destructor
    virtual ~icoMulticomponentVoF();


    // Member Functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Predict the momentum transport
        virtual void momentumTransportPredictor();

        //- Predict thermophysical transport
        virtual void thermophysicalTransportPredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum transport
        virtual void momentumTransportCorrector();

        //- Correct the thermophysical transport
        virtual void thermophysicalTransportCorrector();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const icoMulticomponentVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
