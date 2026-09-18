/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::fv::boussinesqForce

Description
    Must be updated. ---IGNORE---    
    Discrete Tray Radiation Method for collimated radiation flux, derived
    from the implementation given in OpenFOAM-v2512.

    The ray is considered planar and developed to work with a solver
    derived from compressibleVoF.

SourceFiles
    boussinesqForce.C

\*---------------------------------------------------------------------------*/

#ifndef boussinesqForce_H
#define boussinesqForce_H

#include "fvModel.H"

#include "volFields.H"
#include "fluidThermo.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class boussinesqForce Declaration
\*---------------------------------------------------------------------------*/

class boussinesqForce
:
    public fvModel
{
    // Private Data

        // Reference to the phase fraction
        const volScalarField& alpha_;

        // Reference to the temperature field
        const volScalarField& T_;

        //- Reference to the phase thermo
        const fluidThermo& thermo_;

        // Thermal expansion coefficient
        const dimensionedScalar beta_;

        // Reference temperature
        const dimensionedScalar T0_;

        //- Gravitational acceleration
        uniformDimensionedVectorField g_;

public:

    //- Runtime type information
    TypeName("boussinesqForce");


    // Constructors

        //- Construct from explicit source name and mesh
        boussinesqForce
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        boussinesqForce
        (
            const boussinesqForce&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- This does nothing
            virtual void correct();


        // Add explicit and implicit contributions to compressible equation

            //- Add implicit contribution to mixture momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const boussinesqForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
