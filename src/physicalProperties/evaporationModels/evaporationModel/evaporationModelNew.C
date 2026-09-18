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

\*---------------------------------------------------------------------------*/

#include "evaporationModel.H"
#include "noEvaporationModel.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::evaporationModel> Foam::evaporationModel::New
(
    const fvMesh& mesh,
    const word& group
)
{
    const IOdictionary dict
    (
        evaporationModel::findModelDict(mesh, group)
    );

    if (dict.isDict(dictName_))
    {
        const dictionary& modelDict = dict.subDict(dictName_);

        const word modelType(modelDict.lookup("type"));

        Info<< "Selecting evaporation model " << modelType << endl;

        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(modelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown evaporation model " << modelType << nl << nl
                << "Valid evaporation models are : " << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<evaporationModel>(cstrIter()(mesh, group));
    }
    else
    {
        Info<<"There is no " << dictName_ << " dictionary"<<endl;
        Info<<"Selecting default evaporation model none"<<endl;
        return autoPtr<evaporationModel>
        (
            new evaporationModels::none(mesh, group)
        );
    }
}


// ************************************************************************* //