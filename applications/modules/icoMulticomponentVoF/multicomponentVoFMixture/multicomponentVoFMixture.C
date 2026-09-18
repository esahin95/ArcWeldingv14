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

#include "correctContactAngle.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multicomponentVoFMixture, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multicomponentVoFMixture::multicomponentVoFMixture
(
    const fvMesh& mesh
)
:
    compressibleMultiphaseVoFMixture(mesh),

    phases_(phases()),

    missible_(phases().size(), false),

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
    {
        wordList missible(lookup("missible"));
        forAll(phases_, phasei)
        {
            forAll(missible, phasej)
            {
                missible_[phasei] =
                    missible_[phasei] ||
                    (
                        phases_[phasei].name() == missible[phasej]
                    );
            }

            phases_[phasei].vDot().writeOpt() = IOobject::NO_WRITE;
        }
        Info<< "Missible phases: " << missible_ << endl;
    }

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

            sigmaPtrs_.insert
            (
                iter.key(),
                surfaceTensionModel::New
                (
                    dict,
                    mesh
                )
            );
        }
    }

    forAll(phases_, phasei)
    {
        const compressibleVoFphase& alpha1 = phases_[phasei];

        for (label phasej = phasei+1; phasej<phases_.size(); phasej++)
        {
            const compressibleVoFphase& alpha2 = phases_[phasej];

            sigmaPtrTable::const_iterator sigmaPtr =
                sigmaPtrs_.find(interfacePair(alpha1, alpha2));

            if (sigmaPtr == sigmaPtrs_.end() && !missible(phasei, phasej))
            {
                FatalErrorInFunction
                    << "Cannot find interface "
                    << interfacePair(alpha1, alpha2)
                    << " in list of sigma dictionaries"
                    << exit(FatalError);
            }

            sigmaTable::const_iterator D =
                Ds_.find(interfacePair(alpha1, alpha2));

            if (D == Ds_.end() && missible(phasei, phasej))
            {
                FatalErrorInFunction
                    << "Cannot find binary mass diffusion "
                    << interfacePair(alpha1, alpha2)
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
                    debug ?
                    IOobject::AUTO_WRITE :
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimKinematicViscosity, Zero)
            )
        );
    }

    correct();

    updateDm();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

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

            tmp<volScalarField> tsigma
            (
                volScalarField::New("sigma", mesh_, dimSigma_)
            );
            volScalarField& sigma = tsigma.ref();

            sigmaPtrTable::const_iterator sigmaPtr =
                sigmaPtrs_.find(interfacePair(alpha1, alpha2));

            if (sigmaPtr == sigmaPtrs_.end())
            {
                sigma = Zero;
            }
            else
            {
                sigma = sigmaPtr()->sigma();
            }

            const surfaceVectorField gradSigma
            (
                fvc::interpolate(fvc::grad(sigma))
            );

            const surfaceVectorField nHat
            (
                nHatfv(alpha1, alpha2)
            );

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


void Foam::multicomponentVoFMixture::correctThermo()
{
    compressibleMultiphaseVoFMixture::correctThermo();
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
        tkappaEff.ref() +=
            phases_[phasei]
           *(
               phases_[phasei].thermo().kappa()
             + phases_[phasei].thermo().rho()*phases_[phasei].thermo().Cp()*nut
            );
    }

    return tkappaEff;
}

const Foam::volScalarField& Foam::multicomponentVoFMixture::Dm
(
    const label phasei
) const
{
    // volScalarField::New("Dm", rho()*Dm_[phasei])
    return Dm_[phasei];
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


void Foam::multicomponentVoFMixture::updateDm() const
{
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

    dimensionedScalar lBound(sumRhoByWD.dimensions(), 1);

    forAll(phases_, phasei)
    {
        if (!missible_[phasei]) continue;
        const compressibleVoFphase& alpha1 = phases_[phasei];

        sumRhoByWD = Zero;

        forAll(phases_, phasej)
        {
            const compressibleVoFphase& alpha2 = phases_[phasej];
            if (&alpha1 == &alpha2) continue;

            sigmaTable::const_iterator D =
                Ds_.find(interfacePair(alpha1, alpha2));

            dimensionedScalar rDij
            (
                dimless/dimKinematicViscosity,
                D == Ds_.end()? 1e8 : 1.0/max(D(), 1e-8)
            );

            sumRhoByWD += alpha2
                         *alpha2.thermo().rho()
                         /alpha2.thermo().W()
                         *rDij;
        }

        Dm_[phasei] =
            (
                rhoByW - alpha1*alpha1.thermo().rho()/alpha1.thermo().W()
            ) /
            (
                sumRhoByWD + lBound
            );
    }
}

// ************************************************************************* //
