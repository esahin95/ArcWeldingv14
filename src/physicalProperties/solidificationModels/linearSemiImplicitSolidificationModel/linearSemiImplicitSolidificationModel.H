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
    Foam::solidificationModels::linearSemiImplicit

Description
    Solidification with a solid fraction varying linearly between 0 at the
    liquidus and 1 at the solidus temperature; the latent heat is linearised
    in the temperature.

    Within the temperature corrections the solid fraction is updated as
    s = s0 + dsdT*(T - T0), where the slope dsdT and temperature T0 linearise
    the solid fraction about the current state. Outside the mushy zone a
    large slope drives the solid fraction back to 0 or 1. The latent heat
    source is split into a part implicit in T and an explicit remainder. The
    mushy-zone damping and the parameters are those of linearExplicit.

SourceFiles
    linearSemiImplicitSolidificationModel.C

See also
    Foam::solidificationModels::linearExplicit

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
                     Class linearSemiImplicit Declaration
\*---------------------------------------------------------------------------*/

class linearSemiImplicit
:
    public linearExplicit
{
protected:

    // Protected Data

        //- Linearisation temperature
        mutable volScalarField T0_;

        //- Slope of the solid fraction with respect to temperature
        mutable volScalarField dsdT_;


private:

    // Private Member Functions

        //- Set the slope dsdT and temperature T0 linearising the solid
        //  fraction s at temperature T
        void linearise
        (
            const scalar T,
            const scalar s,
            scalar& dsdT,
            scalar& T0
        ) const;


public:

    //- Runtime type information
    TypeName("linearSemiImplicit");


    // Constructors

        //- Construct for the phase given by group
        linearSemiImplicit(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        linearSemiImplicit(const linearSemiImplicit&) = delete;


    //- Destructor
    virtual ~linearSemiImplicit()
    {}


    // Member Functions

        //- Correct the solid fraction and return its maximum change
        virtual scalar correct(const bool relax=true);

        //- Momentum damping as linearExplicit
        using linearExplicit::addSup;

        //- Add the linearised latent heat source to the energy equation
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
