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
    Foam::solidificationModels::linearImplicit

Description
    Solidification with a solid fraction varying linearly between 0 at the
    liquidus and 1 at the solidus temperature; the latent heat is implicit.

    In the mushy zone the latent heat enters the energy equation as the
    apparent heat capacity alpha*rho*L/(Tliq - Tsol), implicit in T. The solid
    fraction is set to its equilibrium value, under-relaxed by relax. The
    mushy-zone damping and the parameters are those of linearExplicit.

SourceFiles
    linearImplicitSolidificationModel.C

See also
    Foam::solidificationModels::linearExplicit

\*---------------------------------------------------------------------------*/

#ifndef linearImplicitSolidificationModel_H
#define linearImplicitSolidificationModel_H

#include "linearExplicitSolidificationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidificationModels
{

/*---------------------------------------------------------------------------*\
                       Class linearImplicit Declaration
\*---------------------------------------------------------------------------*/

class linearImplicit
:
    public linearExplicit
{
protected:

    // Protected Data

        //- Apparent heat capacity of the latent heat
        mutable volScalarField rhoCpLatent_;


public:

    //- Runtime type information
    TypeName("linearImplicit");


    // Constructors

        //- Construct for the phase given by group
        linearImplicit(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        linearImplicit(const linearImplicit&) = delete;


    //- Destructor
    virtual ~linearImplicit()
    {}


    // Member Functions

        //- Correct the solid fraction and return its maximum change
        virtual scalar correct(const bool relax=true);

        //- Momentum damping as linearExplicit
        using linearExplicit::addSup;

        //- Add the latent heat to the energy equation
        virtual void addSup(fvMatrix<scalar>& eqn) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const linearImplicit&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidificationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
