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
    Foam::fv::greyBodyRadiation

Description
    Grey-body radiation from the interface of a VoF phase to surroundings at
    temperature T0, added to the temperature equation.

    The radiative flux eps*sigma*(T0^4 - T^4) is applied per unit interface
    area, |grad(alpha)|, and linearised about the current temperature T*:

        eps*sigma*(T0^4 + 3*T*^4)*|grad(alpha)|
      - 4*eps*sigma*T*^3*|grad(alpha)|*T

    where the last term is implicit in T. The source is written as the field
    \c radiation.

Usage
    \verbatim
    greyBodyRadiation
    {
        type    greyBodyRadiation;
        phase   metal;
        eps     0.4;
        T0      300;
    }
    \endverbatim

    \table
        Property | Description                     | Required | Default
        phase    | Name of the radiating phase     | yes      |
        eps      | Emissivity                      | yes      |
        T0       | Temperature of the surroundings | yes      |
    \endtable

SourceFiles
    greyBodyRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef greyBodyRadiation_H
#define greyBodyRadiation_H

#include "fvModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class greyBodyRadiation Declaration
\*---------------------------------------------------------------------------*/

class greyBodyRadiation
:
    public fvModel
{
    // Private Data

        //- Phase fraction of the radiating phase
        const volScalarField& alpha_;

        //- Emissivity times the Stefan-Boltzmann constant
        const dimensionedScalar eps_;

        //- Temperature of the surroundings
        const dimensionedScalar T0_;

        //- Interface area density, |grad(alpha)|
        volScalarField delta_;

        //- Radiation source, for output
        mutable volScalarField radiation_;


public:

    //- Runtime type information
    TypeName("greyBodyRadiation");


    // Constructors

        //- Construct from explicit source name and mesh
        greyBodyRadiation
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        greyBodyRadiation(const greyBodyRadiation&) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source
            //  term to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- Update the interface area density
            virtual void correct();


        // Add explicit and implicit contributions to compressible equation

            //- Add the linearised radiation to the temperature equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& T,
                fvMatrix<scalar>& eqn
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
        void operator=(const greyBodyRadiation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
