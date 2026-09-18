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
    Foam::evaporationModel

Description

SourceFiles
    evaporationModelI.H
    evaporationModel.C
    evaporationModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef evaporationModel_H
#define evaporationModel_H

#include "physicalProperties.H"
#include "dimensionedTypes.H"
#include "volFields.H"
#include "fvMatrix.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class evaporationModel Declaration
\*---------------------------------------------------------------------------*/

class evaporationModel
:
    public physicalProperties
{
protected:

    //- Evaporation dictionary name
    static const word dictName_;

    // References to database

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the phase thermo
        const fluidThermo& thermo_;


    // References to fields

        //- Fluid fraction
        const volScalarField& alpha_;

        //- Density
        const volScalarField& rho_;

        //- Temperature
        const volScalarField& T_;

        //- Mass flux
        const surfaceScalarField& alphaRhoPhi_;


    //- Mass transfer rate
    volScalarField::Internal mDot_;


public:

    //- Runtime type information
    TypeName("evaporationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            evaporationModel,
            dictionary,
            (
                const fvMesh& mesh,
                const word& group
            ),
            (mesh, group)
        );


    // Constructors

        //- Construct from components
        evaporationModel(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        evaporationModel(const evaporationModel&) = delete;


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<evaporationModel> New
        (
            const fvMesh& mesh,
            const word& group = word::null
        );


    //- Destructor
    virtual ~evaporationModel()
    {}

    // Member Functions

        //- Correct solid fraction
        virtual scalar correct(const bool relax=true) = 0;

        //- Source term for energy equation
        virtual void addSup(fvMatrix<scalar>& eqn) const = 0;

        //- Source term for momentum equation
        virtual void addSup(fvMatrix<vector>& eqn) const = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const evaporationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
