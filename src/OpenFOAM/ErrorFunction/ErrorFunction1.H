/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::Function1s::ErrorFunction

Description
    Solid fraction as a smooth function of temperature, based on the error
    function:

        f(T) = 0.5*(1 - erf(a*(T - Tmid)))

    with Tmid = (Tsol + Tliq)/2 and a = 4/(Tliq - Tsol). f is 1 in the solid,
    0 in the liquid and 0.5 at Tmid; at Tsol and Tliq it is within 0.25% of
    1 and 0 respectively. For Tliq = Tsol it reduces to a step at Tmid.

Usage
    \verbatim
        <name>
        {
            type    errorFunction;
            Tsol    1650;
            Tliq    1700;
        }
    \endverbatim

SourceFiles
    ErrorFunction1.C

\*---------------------------------------------------------------------------*/

#ifndef ErrorFunction1_H
#define ErrorFunction1_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                        Class ErrorFunction Declaration
\*---------------------------------------------------------------------------*/

class ErrorFunction
:
    public FieldFunction1<scalar, ErrorFunction>
{
    // Private Data

        //- Liquidus temperature
        const scalar Tliq_;

        //- Solidus temperature
        const scalar Tsol_;

        //- Mid-point of the solidification range
        const scalar Tmid_;

        //- Steepness, 4/(Tliq - Tsol)
        const scalar a_;


public:

    //- Runtime type information
    TypeName("errorFunction");


    // Constructors

        //- Construct from name and dictionary
        ErrorFunction
        (
            const word& name,
            const unitSets& units,
            const dictionary& dict
        );

        //- Copy constructor
        ErrorFunction(const ErrorFunction& errf);


    //- Destructor
    virtual ~ErrorFunction();


    // Member Functions

        //- Return the solid fraction at temperature x
        virtual scalar value(const scalar x) const;

        //- Integral between two values; not implemented
        virtual scalar integral(const scalar x1, const scalar x2) const;

        //- Write data to dictionary stream
        virtual void write(Ostream& os, const unitSets& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ErrorFunction&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
