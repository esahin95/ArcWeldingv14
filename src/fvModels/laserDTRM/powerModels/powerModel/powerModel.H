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

SourceFiles
    powerModelI.H
    powerModel.C
    powerModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef powerModel_H
#define powerModel_H

#include "IOdictionary.H"
#include "runTimeSelectionTables.H"

#include "Function1.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class powerModel Declaration
\*---------------------------------------------------------------------------*/

class powerModel
{

protected:

    // References

        //- Constant reference to the mesh
        const fvMesh& mesh_;

    // Laser data

        //- Laser position
        const autoPtr<Function1<vector>> pos_;

        //- Laser radius
        const scalar rad_;

        //- Laser normal direction
        const vector normal_;

        //- Laser power
        const scalar Q_;

    // Containers

        //- Positions
        List<vector> positions_;

        //- Powers
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

        //- Construct from components
        powerModel(const dictionary& dict, const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        powerModel(const powerModel&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<powerModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~powerModel()
    {}


    // Member Functions

        // Initialise positions and powers
        virtual void initialise() = 0;

        // Access positions
        const List<vector>& positions() const
        {
            return positions_;
        }

        // Access powers
        const List<scalar>& powers() const
        {
            return powers_;
        }

        // Access normal
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
