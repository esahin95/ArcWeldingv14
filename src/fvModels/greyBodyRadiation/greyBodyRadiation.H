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
    Discrete Tray Radiation Method for collimated radiation flux, derived
    from the implementation given in OpenFOAM-v2512.

    The ray is considered planar and developed to work with a solver
    derived from compressibleVoF.

Usage
    Example usage:
    \verbatim
    greyBodyRadiation
    {
        type    greyBodyRadiation;
        phase   water;
        ***Add more stuff here***
    }
    \endverbatim

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

        // Reference to the liquid phase
        const volScalarField& alpha_;

        // Emissivity times Stefan-Boltzmann constant
        const dimensionedScalar eps_;

        // Reference temperature
        const dimensionedScalar T0_;

        // Numerical delta function
        volScalarField delta_;

        // Radiation field
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
        greyBodyRadiation
        (
            const greyBodyRadiation&
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

            //- Add explicit and implicit contribution to phase energy equation
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
