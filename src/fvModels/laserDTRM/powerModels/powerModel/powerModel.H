/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::powerModel

Description
    Base class for the laser power models of laserDTRM, which distribute the
    laser power over the start positions of the rays.

    The beam is centred on \c position, a Function1 of time, has the given
    \c radius and points in the \c normal direction.

SourceFiles
    powerModel.C
    powerModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef powerModel_H
#define powerModel_H

#include "fvMesh.H"
#include "Function1.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class powerModel Declaration
\*---------------------------------------------------------------------------*/

class powerModel
{
protected:

    // Protected Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Beam centre as a function of time
        const autoPtr<Function1<vector>> pos_;

        //- Beam radius
        const scalar rad_;

        //- Beam direction, normalised
        const vector normal_;

        //- Laser power
        const scalar Q_;

        //- Start positions of the rays
        List<vector> positions_;

        //- Powers of the rays
        List<scalar> powers_;


public:

    //- Runtime type information
    TypeName("powerModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            powerModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        //- Construct from dictionary and mesh
        powerModel(const dictionary& dict, const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        powerModel(const powerModel&) = delete;


    // Selectors

        //- Select from dictionary
        static autoPtr<powerModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~powerModel()
    {}


    // Member Functions

        //- Set the start positions and powers of the rays
        virtual void initialise() = 0;

        //- Return the start positions of the rays
        const List<vector>& positions() const
        {
            return positions_;
        }

        //- Return the powers of the rays
        const List<scalar>& powers() const
        {
            return powers_;
        }

        //- Return the beam direction
        const vector& normal() const
        {
            return normal_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const powerModel&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
