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
    Foam::linearSemiImplicitSolidificationModel

Description

SourceFiles
    linearSemiImplicitSolidificationModel.C

\*---------------------------------------------------------------------------*/

#ifndef linearSemiImplicitSolidificationModel_H
#define linearSemiImplicitSolidificationModel_H

#include "linearExplicitSolidificationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solidificationModels
{


/*---------------------------------------------------------------------------*\
                          Class linearSemiImplicitSolidificationModel Declaration
\*---------------------------------------------------------------------------*/

class linearSemiImplicit
:
    public linearExplicit
{

protected:

    //- Consistent temperature field
    mutable volScalarField T0_;

    //- Solid fraction slope
    mutable volScalarField dsdT_;


public:

    //- Runtime type information
    TypeName("linearSemiImplicit");

    // Constructors

        //- Construct from components
        linearSemiImplicit(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        linearSemiImplicit(const linearSemiImplicit&) = delete;


    //- Destructor
    virtual ~linearSemiImplicit()
    {}


    // Member functions

        //- Correct solid fraction
        virtual scalar correct(const bool relax=true);

        //- Source term for energy equation
        virtual void addSup(fvMatrix<scalar>& eqn) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const linearSemiImplicit&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidificationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
