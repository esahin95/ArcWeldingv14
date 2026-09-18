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

\*---------------------------------------------------------------------------*/

#ifndef ErrorFunction1_H
#define ErrorFunction1_H

#include "Function1.H"

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
    
    // Private data
    
        const scalar Tliq_;
        const scalar Tsol_;
        const scalar Tmid_;
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

    // Destructor
    
        virtual ~ErrorFunction();

    
    // Member functions

        //- Return error function value as a function of scalar x
        virtual scalar value(const scalar x) const;

        //- Integrate between two scalar fields
        virtual scalar integral(const scalar x1, const scalar x2) const;
    
        //- Write in dictionary format
        virtual void write(Ostream& os, const unitSets& units) const;

    // Member operators

        //- Disallow default bitwise assignment
        void operator=(const ErrorFunction&) = delete;

};

} // End namespace Function1s
} // End namespace Foam


#endif