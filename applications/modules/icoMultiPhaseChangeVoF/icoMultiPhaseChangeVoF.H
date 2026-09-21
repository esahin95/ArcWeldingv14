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
    Foam::solvers::icoMultiPhaseChangeVoF

Description
    Solver module for multiple incompressible, non-isothermal fluids using a
    VOF (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing. Phases can be declared miscible: between miscible phases there
    is no surface tension and no interface compression, but mass diffusion.

    Inherits from icoMulticomponentVoF and adds the runtime selectable
    solidification and evaporation phase change models of every phase. The
    temperature equation is iterated with the phase change models until they
    have converged.

    The numerical controls are given in the PIMPLE dictionary:
    \table
        Property    | Description                           | Default
        nThermoCorr | Maximum number of temperature corrections | 20
        thermoTol   | Tolerance on the change of the phase change models \\
                                                            | 1e-6
    \endtable

SourceFiles
    icoMultiPhaseChangeVoF.C
    thermophysicalPredictor.C
    momentumPredictor.C

See also
    Foam::solvers::icoMulticomponentVoF
    Foam::solvers::icoPhaseChangeVoF

\*---------------------------------------------------------------------------*/

#ifndef icoMultiPhaseChangeVoF_H
#define icoMultiPhaseChangeVoF_H

#include "icoMulticomponentVoF.H"
#include "solidificationModel.H"
#include "evaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                   Class icoMultiPhaseChangeVoF Declaration
\*---------------------------------------------------------------------------*/

class icoMultiPhaseChangeVoF
:
    public icoMulticomponentVoF
{
    // Phase change models

        //- Solidification models of the phases
        PtrList<solidificationModel> solidificationModels_;

        //- Evaporation models of the phases
        PtrList<evaporationModel> evaporationModels_;


    // Numerical Controls

        //- Maximum number of temperature corrections per time step
        const label nThermoCorr_;

        //- Tolerance on the change of the phase change models
        const scalar thermoTol_;


    // Private Member Functions

        //- Add the phase change sources to the equation
        template<class Type>
        void addSup(fvMatrix<Type>& eqn) const;

        //- Correct the phase change models and return their maximum change
        scalar correctPhaseChange();


public:

    //- Runtime type information
    TypeName("icoMultiPhaseChangeVoF");


    // Constructors

        //- Construct from region mesh
        icoMultiPhaseChangeVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        icoMultiPhaseChangeVoF(const icoMultiPhaseChangeVoF&) = delete;


    //- Destructor
    virtual ~icoMultiPhaseChangeVoF();


    // Member Functions

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the temperature equation, iterating with the
        //  phase change models until converged
        virtual void thermophysicalPredictor();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const icoMultiPhaseChangeVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "icoMultiPhaseChangeVoFTemplates.C"
#endif

#endif

// ************************************************************************* //
