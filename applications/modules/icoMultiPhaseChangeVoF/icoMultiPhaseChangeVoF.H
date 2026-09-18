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
    Solver module for the solution of multiple incompressible, non-isothermal
    immiscible fluids using a VOF (volume of fluid) phase-fraction based
    interface capturing approach, with optional mesh motion and mesh topology
    changes including adaptive re-meshing. Selection of phases can be set to
    be missible with zero surface tension and interface compression, but
    non-zero diffusivity.

    Inherits from icoMulticomponentVoF and adds solidification and evaporation
    phase change models for first given phase. Models are runtime selectable.
    Energy equation is iterated until given tolerance is reached.

SourceFiles
    icoMultiPhaseChangeVoF.C

See also
    Foam::solvers::icoMulticomponentVoF

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
protected:

    // Phase change models

        // Solidification model pointers
        PtrList<solidificationModel> solModels_;

        // Evaporation model pointers
        PtrList<evaporationModel> evaModels_;

    // Numerical Controls

        // Thermo corrections
        const label nThermoCorr_;

        // Thermo tolerance
        const scalar thermoTol_;


    // Private member functions

        // Add source terms
        template<class Type>
        void addSup(fvMatrix<Type>& eqn) const;

        // Correct phase change models
        scalar correctPhaseChange();

public:

    //- Runtime type information
    TypeName("icoMultiPhaseChangeVoF");


    // Constructors

        //- Construct from region mesh
        icoMultiPhaseChangeVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        icoMultiPhaseChangeVoF
        (
            const icoMultiPhaseChangeVoF&
        ) = delete;


    //- Destructor
    virtual ~icoMultiPhaseChangeVoF();


    // Member Functions

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation
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
