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
    Compressible multiphase mixture for interface-capturing simulations.

    Derived from viscosity so that it can be used in conjunction with the
    compressible turbulence models.

    Surface tension and contact-angle is handled for the interface between each
    phase-pair. Surface tension and interface compression are ignored between
    missible phases, which have diffusion instead.

SourceFiles
    multicomponentVoFMixture.C

See also
    Foam::compressibleVoFMixture

\*---------------------------------------------------------------------------*/

#ifndef multicomponentVoFMixture_H
#define multicomponentVoFMixture_H

#include "compressibleMultiphaseVoFMixture.H"

#include "surfaceTensionModel.H"
#include "Function2.H"

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
    // Private Data

        //- Phases
        UPtrListDictionary<compressibleVoFphase>& phases_;

        //- Components
        List<bool> missible_;

        //- Heat capacity
        volScalarField rhoCp_;

        //- surface tension models
        typedef HashTable
            <
                dictionary,
                interfacePair,
                interfacePair::hash
            > dictTable;
        typedef HashTable
            <
                autoPtr<surfaceTensionModel>,
                interfacePair,
                interfacePair::hash
            > sigmaPtrTable;
        sigmaPtrTable sigmaPtrs_;

        //- Binary mass diffusion coefficient functions
        sigmaTable Ds_;

        //- Mass diffusion coefficients in the mixture
        mutable PtrList<volScalarField> Dm_;


public:

    TypeName("multicomponentVoFMixture");


    // Constructors

        //- Construct from fvMesh
        multicomponentVoFMixture(const fvMesh& mesh);


    //- Destructor
    virtual ~multicomponentVoFMixture()
    {}


    // Member Functions

        //- Surface tension force
        tmp<surfaceScalarField> surfaceTensionForce
        (
            const volVectorField& U
        ) const;

        //- Access heat capacity
        const volScalarField& rhoCp() const
        {
            return rhoCp_;
        }

        //- Correct the thermodynamics of each phase
        virtual void correctThermo();

        //- Update properties
        virtual void correct();

        //- Return effective heat conductivity for the mixture
        virtual tmp<volScalarField> kappaEff
        (
            const volScalarField& nut
        ) const;

        //- Phases missible?
        bool missible(label phasei, label phasej) const
        {
            return missible_[phasei] && missible_[phasej];
        }

        //- Mixture mass diffusion for phase
        //tmp<volScalarField> DEff(const label phasei) const;
        const volScalarField& Dm(const label phasei) const;

        //- Mass diffusion flux for phase
        tmp<surfaceScalarField> j(const label phasei) const;

        //- Update mixture diffusion coefficients
        void updateDm() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
