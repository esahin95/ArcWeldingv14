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

\*---------------------------------------------------------------------------*/

#include "multicomponentVoFMixture.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multicomponentVoFMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multicomponentVoFMixture::multicomponentVoFMixture(const fvMesh& mesh)
:
    compressibleMultiphaseVoFMixture(mesh),

    phases_(phases()),

    miscible_(phases().size(), false),

    rhoCp_
    (
        IOobject
        (
            "rhoCp",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar("rhoCp", dimEnergy/dimVolume/dimTemperature, 0)
    ),

    Ds_(lookup("Ds")),

    Dm_(phases().size())
{
    const wordList miscibleNames
    (
        lookupBackwardsCompatible<wordList>({"miscible", "missible"})
    );

    forAll(phases_, phasei)
    {
        miscible_[phasei] =
            findIndex(miscibleNames, phases_[phasei].name()) != -1;

        phases_[phasei].vDot().writeOpt() = IOobject::NO_WRITE;
    }

    Info<< "Miscible phases: " << miscible_ << endl;

    if (found("sigmaDicts"))
    {
        const dictTable sigmaDicts(lookup("sigmaDicts"));

        forAllConstIter(dictTable, sigmaDicts, iter)
        {
            sigmaPtrs_.insert
            (
                iter.key(),
                surfaceTensionModel::New(iter(), mesh)
            );
        }
    }
    else
    {
        forAllConstIter(sigmaTable, sigmas_, iter)
        {
            dictionary dict;
            dict.add("sigma", iter());

            sigmaPtrs_.insert(iter.key(), surfaceTensionModel::New(dict, mesh));
        }
    }

    // Check that every immiscible pair has a surface tension and every
    // miscible pair a diffusion coefficient
    forAll(phases_, phasei)
    {
        const compressibleVoFphase& alpha1 = phases_[phasei];

        for (label phasej = phasei+1; phasej<phases_.size(); phasej++)
        {
            const compressibleVoFphase& alpha2 = phases_[phasej];
            const interfacePair pair(alpha1, alpha2);

            if (!sigmaPtrs_.found(pair) && !miscible(phasei, phasej))
            {
                FatalErrorInFunction
                    << "Cannot find interface " << pair
                    << " in list of sigma dictionaries"
                    << exit(FatalError);
            }

            if (!Ds_.found(pair) && miscible(phasei, phasej))
            {
                FatalErrorInFunction
                    << "Cannot find binary mass diffusion " << pair
                    << " in list of interfaces"
                    << exit(FatalError);
            }
        }
    }

    forAll(phases_, phasei)
    {
        Dm_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("Dm", phases_[phasei].name()),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    debug ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimKinematicViscosity, Zero)
            )
        );
    }

    correct();

    updateDm();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::multicomponentVoFMixture::surfaceTensionForce
(
    const volVectorField& U
) const
{
    tmp<surfaceScalarField> tstf
    (
        surfaceScalarField::New
        (
            "surfaceTensionForce",
            mesh_,
            dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), 0)
        )
    );

    surfaceScalarField& stf = tstf.ref();

    forAll(phases_, phasei)
    {
        const compressibleVoFphase& alpha1 = phases_[phasei];

        for (label phasej = phasei+1; phasej<phases_.size(); phasej++)
        {
            const compressibleVoFphase& alpha2 = phases_[phasej];

            sigmaPtrTable::const_iterator sigmaPtr =
                sigmaPtrs_.find(interfacePair(alpha1, alpha2));

            // No surface tension between miscible phases
            if (sigmaPtr == sigmaPtrs_.end())
            {
                continue;
            }

            tmp<volScalarField> tsigma
            (
                volScalarField::New("sigma", mesh_, dimSigma_)
            );
            volScalarField& sigma = tsigma.ref();
            sigma = sigmaPtr()->sigma();

            const surfaceVectorField gradSigma
            (
                fvc::interpolate(fvc::grad(sigma))
            );

            const surfaceVectorField nHat(nHatfv(alpha1, alpha2));

            // Normal force from the curvature, and tangential (Marangoni)
            // force from the surface tension gradient along the interface
            stf += fvc::interpolate(sigma*K(alpha1, alpha2, U))*
                (
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                ) +
                (
                    (
                        fvc::interpolate
                        (
                            mag
                            (
                                alpha2*fvc::grad(alpha1)
                              - alpha1*fvc::grad(alpha2)
                            )
                        ) * (gradSigma - (gradSigma&nHat)*nHat)
                    ) & (mesh_.Sf() / mesh_.magSf())
                );
        }
    }

    return tstf;
}


void Foam::multicomponentVoFMixture::correct()
{
    compressibleMultiphaseVoFMixture::correct();

    rhoCp_ = Zero;
    forAll(phases_, phasei)
    {
        rhoCp_ += phases_[phasei]
                 *phases_[phasei].thermo().rho()
                 *phases_[phasei].thermo().Cp();
    }
}


Foam::tmp<Foam::volScalarField> Foam::multicomponentVoFMixture::kappaEff
(
    const volScalarField& nut
) const
{
    tmp<volScalarField> tkappaEff
    (
        phases_[0]
       *(
            phases_[0].thermo().kappa()
          + phases_[0].thermo().rho()*phases_[0].thermo().Cp()*nut
        )
    );

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        const rhoFluidThermo& thermo = phases_[phasei].thermo();

        tkappaEff.ref() +=
            phases_[phasei]*(thermo.kappa() + thermo.rho()*thermo.Cp()*nut);
    }

    return tkappaEff;
}


Foam::tmp<Foam::surfaceScalarField> Foam::multicomponentVoFMixture::j
(
    const label phasei
) const
{
    return surfaceScalarField::New
    (
        "j",
        -fvc::interpolate(Dm_[phasei])
        *mesh_.magSf()
        *fvc::snGrad(phases_[phasei])
    );
}


void Foam::multicomponentVoFMixture::updateDm()
{
    // Total molar concentration
    volScalarField rhoByW
    (
        volScalarField::New
        (
            "rhoByW",
            mesh_,
            dimensionedScalar(dimMoles/dimVolume, Zero)
        )
    );

    forAll(phases_, phasei)
    {
        rhoByW += phases_[phasei]
                 *phases_[phasei].thermo().rho()
                 /phases_[phasei].thermo().W();
    }

    volScalarField sumRhoByWD
    (
        volScalarField::New
        (
            "sumRhoByWD",
            mesh_,
            dimMoles/dimVolume/dimKinematicViscosity
        )
    );

    // Regularisation of the denominator where a phase is on its own
    const dimensionedScalar lBound(sumRhoByWD.dimensions(), 1);

    forAll(phases_, phasei)
    {
        if (!miscible_[phasei]) continue;

        const compressibleVoFphase& alpha1 = phases_[phasei];

        sumRhoByWD = Zero;

        forAll(phases_, phasej)
        {
            const compressibleVoFphase& alpha2 = phases_[phasej];

            if (&alpha1 == &alpha2) continue;

            sigmaTable::const_iterator D =
                Ds_.find(interfacePair(alpha1, alpha2));

            const dimensionedScalar rDij
            (
                dimless/dimKinematicViscosity,
                D == Ds_.end() ? 1e8 : 1.0/max(D(), 1e-8)
            );

            sumRhoByWD += alpha2
                         *alpha2.thermo().rho()
                         /alpha2.thermo().W()
                         *rDij;
        }

        Dm_[phasei] =
            (rhoByW - alpha1*alpha1.thermo().rho()/alpha1.thermo().W())
           /(sumRhoByWD + lBound);
    }
}


// ************************************************************************* //
