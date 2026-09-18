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
    Foam::solidificationModel

Description
    Base class for solidification models of a VoF phase.

    The model is selected from the optional \c solidification sub-dictionary
    of the phase physicalProperties; without it no solidification is
    modelled (\c none). A model maintains the solid fraction of the phase,
    \c solid.<phase>, and provides the corresponding source terms for the
    energy and momentum equations.

Usage
    In constant/physicalProperties.<phase>:
    \verbatim
    solidification
    {
        type    linearImplicit;
        L       2.7e5;
        Tsol    1650;
        Tliq    1700;
    }
    \endverbatim

SourceFiles
    solidificationModel.C
    solidificationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef solidificationModel_H
#define solidificationModel_H

#include "physicalProperties.H"
#include "fluidThermo.H"
#include "volFields.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class solidificationModel Declaration
\*---------------------------------------------------------------------------*/

class solidificationModel
:
    public physicalProperties
{
protected:

    // Protected Data

        //- Name of the model sub-dictionary
        static const word dictName_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the phase thermo
        const fluidThermo& thermo_;

        //- Phase fraction
        const volScalarField& alpha_;

        //- Temperature
        const volScalarField& T_;

        //- Phase mass flux
        const surfaceScalarField& alphaRhoPhi_;

        //- Solid fraction of the phase
        volScalarField sf_;


public:

    //- Runtime type information
    TypeName("solidificationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidificationModel,
            dictionary,
            (
                const fvMesh& mesh,
                const word& group
            ),
            (mesh, group)
        );


    // Constructors

        //- Construct for the phase given by group
        solidificationModel(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        solidificationModel(const solidificationModel&) = delete;


    // Selectors

        //- Select the solidification model of the phase given by group
        static autoPtr<solidificationModel> New
        (
            const fvMesh& mesh,
            const word& group = word::null
        );


    //- Destructor
    virtual ~solidificationModel()
    {}


    // Member Functions

        //- Correct the solid fraction and return its maximum change
        virtual scalar correct(const bool relax=true) = 0;

        //- Add the latent heat source to the energy equation
        virtual void addSup(fvMatrix<scalar>& eqn) const = 0;

        //- Add the solid damping source to the momentum equation
        virtual void addSup(fvMatrix<vector>& eqn) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidificationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
