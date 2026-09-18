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
    Foam::reflectionModel

Description

SourceFiles
    reflectionModelI.H
    reflectionModel.C
    reflectionModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef reflectionModel_H
#define reflectionModel_H

#include "IOdictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class reflectionModel Declaration
\*---------------------------------------------------------------------------*/

class reflectionModel
{
    // Private Data

        //- Constant reference to the mesh
        const fvMesh& mesh_;

public:

    //- Runtime type information
    TypeName("reflectionModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            reflectionModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        //- Construct from components
        reflectionModel(const dictionary& dict, const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        reflectionModel(const reflectionModel&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<reflectionModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~reflectionModel()
    {}


    // Member Functions

        // Reflected direction
        virtual vector R
        (
            const vector& incident,
            const vector& n
        ) const;

        // Reflectivity
        virtual scalar rho(const scalar cosTheta) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reflectionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
