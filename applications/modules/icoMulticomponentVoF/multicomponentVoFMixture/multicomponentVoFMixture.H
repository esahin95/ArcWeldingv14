/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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
    Foam::multicomponentVoFMixture

Description
    Incompressible multiphase mixture for interface-capturing simulations in
    which some of the phases may be miscible.

    Surface tension, including the Marangoni force of a temperature dependent
    surface tension, and contact angles are handled for the interface between
    each pair of immiscible phases. Between miscible phases there is no
    surface tension and no interface compression, but mass diffusion with the
    flux j_i = -D_i*grad(alpha_i). The mixture-averaged diffusion coefficient
    D_i of phase i follows from the binary coefficients D_ij and the molar
    concentrations c_j = alpha_j*rho_j/W_j:

        D_i = (c - c_i)/(sum_{j != i}(c_j/D_ij) + 1)

    where c is the total molar concentration, and the 1 [mol s/m^5] avoids
    the division by zero where phase i is on its own. Binary coefficients are
    limited to at least 1e-8 m^2/s, which is also used for the pairs without
    one. The diffusion coefficients are written as Dm.<phase> in debug mode.

Usage
    Entries of constant/phaseProperties in addition to those of
    compressibleMultiphaseVoFMixture:
    \table
        Property   | Description                           | Required
        miscible   | Names of the miscible phases          | yes
        Ds         | Binary diffusion coefficients of the miscible phase \\
                     pairs                                 | yes
        sigmaDicts | Surface tension models of the phase pairs, replacing \\
                     the constant values of sigmas         | no
    \endtable

    Example:
    \verbatim
    phases      (water oil air);

    miscible    (water oil);

    sigmas
    (
        (air water)     0.07
        (air oil)       0.07
        (water oil)     0
    );

    Ds
    (
        (water oil)     1e-5
    );
    \endverbatim

    The keyword missible is accepted in place of miscible.

SourceFiles
    multicomponentVoFMixture.C

See also
    Foam::compressibleMultiphaseVoFMixture

\*---------------------------------------------------------------------------*/

#ifndef multicomponentVoFMixture_H
#define multicomponentVoFMixture_H

#include "compressibleMultiphaseVoFMixture.H"
#include "surfaceTensionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class multicomponentVoFMixture Declaration
\*---------------------------------------------------------------------------*/

class multicomponentVoFMixture
:
    public compressibleMultiphaseVoFMixture
{
    // Private Typedefs

        //- Table of surface tension models of the phase pairs
        typedef HashTable
        <
            autoPtr<surfaceTensionModel>,
            interfacePair,
            interfacePair::hash
        > sigmaPtrTable;

        //- Table of surface tension model dictionaries of the phase pairs
        typedef HashTable
        <
            dictionary,
            interfacePair,
            interfacePair::hash
        > dictTable;


    // Private Data

        //- Phases
        UPtrListDictionary<compressibleVoFphase>& phases_;

        //- Is each phase miscible?
        List<bool> miscible_;

        //- Volumetric heat capacity
        volScalarField rhoCp_;

        //- Surface tension models of the immiscible phase pairs
        sigmaPtrTable sigmaPtrs_;

        //- Binary mass diffusion coefficients of the miscible phase pairs
        sigmaTable Ds_;

        //- Mixture-averaged mass diffusion coefficients of the phases
        PtrList<volScalarField> Dm_;


public:

    //- Runtime type information
    TypeName("multicomponentVoFMixture");


    // Constructors

        //- Construct from fvMesh
        multicomponentVoFMixture(const fvMesh& mesh);


    //- Destructor
    virtual ~multicomponentVoFMixture()
    {}


    // Member Functions

        //- Return the surface tension force of the immiscible interfaces
        tmp<surfaceScalarField> surfaceTensionForce
        (
            const volVectorField& U
        ) const;

        //- Return the volumetric heat capacity
        const volScalarField& rhoCp() const
        {
            return rhoCp_;
        }

        //- Update properties
        virtual void correct();

        //- Return the effective thermal conductivity of the mixture
        virtual tmp<volScalarField> kappaEff
        (
            const volScalarField& nut
        ) const;

        //- Are the given phases both miscible?
        bool miscible(const label phasei, const label phasej) const
        {
            return miscible_[phasei] && miscible_[phasej];
        }

        //- Return the mass diffusion coefficient of the given phase
        const volScalarField& Dm(const label phasei) const
        {
            return Dm_[phasei];
        }

        //- Return the mass diffusion flux of the given phase
        tmp<surfaceScalarField> j(const label phasei) const;

        //- Update the mass diffusion coefficients
        void updateDm();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
