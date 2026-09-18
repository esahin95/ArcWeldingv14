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
    Foam::solidificationModels::linearExplicit

Description
    Solidification with a solid fraction varying linearly between 0 at the
    liquidus and 1 at the solidus temperature; the latent heat is explicit.

    The latent heat release L*D(alpha*rho*s)/Dt is added explicitly to the
    energy equation. Within the temperature corrections the solid fraction
    is moved towards its equilibrium value s_eq by relax*St*(s_eq - s), with
    the Stefan number St = Cp*(Tliq - Tsol)/L. The flow is damped in the
    mushy zone by the Carman-Kozeny sink Cu*a^2/((1 - a)^3 + q) on U, where
    a = alpha*s is written as \c alpha.solid.<phase>.

Usage
    \table
        Property | Description                       | Required | Default
        L        | Latent heat of fusion             | yes      |
        Tliq     | Liquidus temperature              | yes      |
        Tsol     | Solidus temperature               | yes      |
        Cu       | Mushy-zone damping constant       | no       | 1e5
        q        | Damping regularisation constant   | no       | 1e-3
        relax    | Solid fraction relaxation factor  | no       | 1
    \endtable

SourceFiles
    linearExplicitSolidificationModel.C

\*---------------------------------------------------------------------------*/

#ifndef linearExplicitSolidificationModel_H
#define linearExplicitSolidificationModel_H

#include "solidificationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidificationModels
{

/*---------------------------------------------------------------------------*\
                       Class linearExplicit Declaration
\*---------------------------------------------------------------------------*/

class linearExplicit
:
    public solidificationModel
{
protected:

    // Protected Data

        //- Model coefficients dictionary
        const dictionary& dict_;

        //- Latent heat of fusion
        const dimensionedScalar Lm_;

        //- Liquidus temperature
        const dimensionedScalar Tliq_;

        //- Solidus temperature
        const dimensionedScalar Tsol_;

        //- Mushy-zone damping constant
        const dimensionedScalar Cu_;

        //- Damping regularisation constant
        const scalar q_;

        //- Solid fraction relaxation factor
        const scalar relax_;

        //- Solid volume fraction, alpha*s
        volScalarField alphaSolid_;


    // Protected Member Functions

        //- Return the equilibrium solid fraction at the current temperature
        tmp<volScalarField> sfEquilibrium() const;


public:

    //- Runtime type information
    TypeName("linearExplicit");


    // Constructors

        //- Construct for the phase given by group
        linearExplicit(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        linearExplicit(const linearExplicit&) = delete;


    //- Destructor
    virtual ~linearExplicit()
    {}


    // Member Functions

        //- Correct the solid fraction and return its maximum change
        virtual scalar correct(const bool relax=true);

        //- Add the latent heat source to the energy equation
        virtual void addSup(fvMatrix<scalar>& eqn) const;

        //- Add the mushy-zone damping to the momentum equation
        virtual void addSup(fvMatrix<vector>& eqn) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const linearExplicit&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidificationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
