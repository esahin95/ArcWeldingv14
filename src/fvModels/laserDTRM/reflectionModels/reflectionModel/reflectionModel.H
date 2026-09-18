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
    Base class for the reflection models of laserDTRM, providing the
    reflectivity as a function of the cosine of the angle of incidence and
    the direction of the specularly reflected ray.

SourceFiles
    reflectionModel.C
    reflectionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef reflectionModel_H
#define reflectionModel_H

#include "fvMesh.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class reflectionModel Declaration
\*---------------------------------------------------------------------------*/

class reflectionModel
{
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

        //- Construct from dictionary and mesh
        reflectionModel(const dictionary& dict, const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        reflectionModel(const reflectionModel&) = delete;


    // Selectors

        //- Select from dictionary
        static autoPtr<reflectionModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~reflectionModel()
    {}


    // Member Functions

        //- Return the direction of the specularly reflected ray for the
        //  incident direction and the unit interface normal n
        virtual vector R(const vector& incident, const vector& n) const;

        //- Return the reflectivity for the cosine of the angle of incidence
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
