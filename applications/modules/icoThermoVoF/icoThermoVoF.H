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

Class
    Foam::solvers::icoThermoVoF

Description
    Solver module for 2 incompressible non-isothermal immiscible fluids using
    a VOF (volume of fluid) phase fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

    The momentum and other fluid properties are of the "mixture_" and a single
    momentum equation is solved.

    Either mixture_ or two-phase transport modelling may be selected.  In the
    mixture_ approach a single laminar, RAS or LES model is selected to model the
    momentum stress.  In the Euler-Euler two-phase approach separate laminar,
    RAS or LES selected models are selected for each of the phases.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources like laser heat source, grey
    body radiation and Boussinesq force.

    Energy equation is build based on icoReactingMultiphaseInterFoam. Alpha
    predictor and pressure corrector are taken from twoPhaseVoFSolver and are
    identical to incompressibleVoF.

SourceFiles
    icoThermoVoF.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::twoPhaseVoFSolver

\*---------------------------------------------------------------------------*/

#ifndef icoThermoVoF_H
#define icoThermoVoF_H

#include "twoPhaseVoFSolver.H"
#include "compressibleTwoPhaseVoFMixture.H"
#include "compressibleInterPhaseTransportModel.H"
#include "compressibleInterPhaseThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class icoThermoVoF Declaration
\*---------------------------------------------------------------------------*/

class icoThermoVoF
:
    public twoPhaseVoFSolver
{

protected:

    // Phase properties

        //- The compressible two-phase mixture_
        compressibleTwoPhaseVoFMixture& mixture_;


    // Thermophysical properties

        //- Reference to mixture pressure field
        volScalarField& p;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference_;


    // Kinematic properties

        // Phase-1 mass-flux
        surfaceScalarField alphaRhoPhi1;

        // Phase-2 mass-flux
        surfaceScalarField alphaRhoPhi2;

    //- Thermophysical properties

        //- Volumetric heat capacity
        volScalarField rhoCp;

        //- Specific heat capacity
        volScalarField Cp;

        //- Heat capacity density flux
        surfaceScalarField rhoPhiCp;


    // Momentum transport

        //- Momentum transport model
        compressibleInterPhaseTransportModel momentumTransport;


    // Thermophysical transport

        //- Thermophysical transport model
        compressibleInterPhaseThermophysicalTransportModel
            thermophysicalTransport;


    // Protected Member Functions

        //- Return the pressure reference
        virtual const Foam::pressureReference& pressureReference() const
        {
            return pressureReference_;
        }

        //- The flow is incompressible
        virtual bool incompressible() const
        {
            return mixture_.incompressible();;
        }

        //- Is the flow divergent?
        //  i.e. includes phase-fraction sources
        virtual bool divergent() const
        {
            return
               !incompressible()
             || fvModels().addsSupToField(alpha1.name())
             || fvModels().addsSupToField(alpha2.name());
        }

        //- Return the mixture_ compressibility/density
        //  Not required for incompressible fluids
        virtual tmp<volScalarField> psiByRho() const
        {
            NotImplemented;
            return mixture_.psiByRho();
        }

        //- Calculate the alpha equation sources
        virtual void alphaSuSp
        (
            tmp<volScalarField::Internal>& Su,
            tmp<volScalarField::Internal>& Sp
        );

        //- Return the momentum equation stress term
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U)
        {
            return momentumTransport.divDevTau(U);
        }


public:

    //- Runtime type information
    TypeName("icoThermoVoF");


    // Constructors

        //- Construct from region mesh
        icoThermoVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        icoThermoVoF(const icoThermoVoF&) = delete;


    //- Destructor
    virtual ~icoThermoVoF();


    // Member Functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Predict the momentum transport
        virtual void momentumTransportPredictor();

        //- Predict thermophysical transport
        virtual void thermophysicalTransportPredictor();

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
        void operator=(const icoThermoVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
