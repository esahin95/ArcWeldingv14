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
    Solver module for multiple incompressible, non-isothermal fluids using a
    VOF (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing. Phases can be declared miscible: between miscible phases there
    is no surface tension and no interface compression, but mass diffusion.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved. A mixture approach for momentum transport is
    provided in which a single laminar, RAS or LES model is selected to model
    the momentum stress.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    The temperature equation is based on icoReactingMultiphaseInterFoam. The
    mass diffusion between miscible phases uses mixture-averaged diffusion
    coefficients and Fick's law, and is added explicitly to the phase
    fraction fluxes. The alpha predictor and pressure corrector are taken
    from incompressibleMultiphaseVoF.

    The time step is limited by the Courant numbers of the base solver and by
    two diffusion numbers, given in the controlDict:
    \table
        Property | Description                                 | Required
        maxEnCo  | Maximum thermal diffusion number            | yes
        maxDiCo  | Maximum mass diffusion number               | yes
    \endtable

    The continuity error is written as the field \c contErr.

SourceFiles
    icoMulticomponentVoF.C
    alphaPredictor.C
    thermophysicalPredictor.C
    pressureCorrector.C

See also
    Foam::solvers::multiphaseVoFSolver
    Foam::multicomponentVoFMixture

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

        //- The multicomponent mixture
        multicomponentVoFMixture& mixture;

        //- Reference to the phases
        UPtrListDictionary<compressibleVoFphase>& phases;


    // Thermophysical properties

        //- Reference to the mixture static pressure field
        volScalarField& p;

        //- Volumetric heat capacity of the mixture
        const volScalarField& rhoCp;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference_;


    // Kinematic properties

        //- Phase mass fluxes, also looked up by the phase change models
        PtrList<surfaceScalarField> alphaRhoPhis_;


    // Momentum transport

        //- Momentum transport model pointer
        autoPtr<compressible::momentumTransportModel> momentumTransport_;

        //- Momentum transport model
        compressible::momentumTransportModel& momentumTransport;


    // Diagnostics

        //- Continuity error, written for post-processing
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
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U)
        {
            return momentumTransport.divDevTau(U);
        }

        //- Solve for the phase-fractions
        virtual void alphaPredictor();

        //- Return the interface surface tension force for the momentum
        //  equation. Overridden as the multicomponent mixture's force, which
        //  omits the miscible interfaces, hides that of multiphaseVoFMixture.
        virtual tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Return the temperature equation without source terms
        tmp<fvScalarMatrix> temperatureEqn();


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
        icoMulticomponentVoF(const icoMulticomponentVoF&) = delete;


    //- Destructor
    virtual ~icoMulticomponentVoF();


    // Member Functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Predict the momentum transport
        virtual void momentumTransportPredictor();

        //- Predict thermophysical transport; does nothing
        virtual void thermophysicalTransportPredictor();

        //- Construct and solve the temperature equation
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum transport
        virtual void momentumTransportCorrector();

        //- Correct the thermophysical transport; does nothing
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
