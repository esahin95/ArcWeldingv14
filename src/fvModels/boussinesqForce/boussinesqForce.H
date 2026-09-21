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
    Boussinesq buoyancy force on a constant-density VoF phase, added to the
    mixture momentum equation.

    The density of the phase is taken to vary with temperature as
    rho*(1 - beta*(T - T0)). The momentum equation already contains the
    gravity force on the constant density rho, so only the variation is
    added:

        -alpha*rho*beta*(T - T0)*g

    where alpha is the phase fraction, beta the thermal expansion coefficient
    and T0 the reference temperature. The phase must have a constant-density
    equation of state, and g is read from constant/g.

Usage
    \verbatim
    boussinesqForce
    {
        type    boussinesqForce;
        phase   metal;
        beta    5e-5;
        T0      1700;
    }
    \endverbatim

    \table
        Property | Description                   | Required | Default
        phase    | Name of the buoyant phase     | yes      |
        beta     | Thermal expansion coefficient | yes      |
        T0       | Reference temperature         | yes      |
    \endtable

SourceFiles
    boussinesqForce.C

\*---------------------------------------------------------------------------*/

#ifndef boussinesqForce_H
#define boussinesqForce_H

#include "fvModel.H"
#include "fluidThermo.H"
#include "volFields.H"
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

        //- Phase fraction
        const volScalarField& alpha_;

        //- Temperature
        const volScalarField& T_;

        //- Phase thermo
        const fluidThermo& thermo_;

        //- Thermal expansion coefficient
        const dimensionedScalar beta_;

        //- Reference temperature
        const dimensionedScalar T0_;

        //- Gravitational acceleration
        const uniformDimensionedVectorField g_;


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
        boussinesqForce(const boussinesqForce&) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source
            //  term to the transport equation
            virtual wordList addSupFields() const;


        // Add explicit and implicit contributions to compressible equation

            //- Add the buoyancy force to the mixture momentum equation
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
