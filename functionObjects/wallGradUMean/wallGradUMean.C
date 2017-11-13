/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "wallGradUMean.H"
#include "turbulenceModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallGradUMean, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        wallGradUMean,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::wallGradUMean::writeFileHeader(const label i)
{
    writeHeader(file(), "y+ ()");

    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "average");
    file() << endl;
}


void Foam::functionObjects::wallGradUMean::calcWallGradUMean
(
    const turbulenceModel& turbModel,
    volVectorField& wallGradUMean
)
{
    volScalarField::Boundary d = nearWallDist(mesh_).y();

    const volScalarField::Boundary nutBf =
        turbModel.nut()().boundaryField();

    const volScalarField::Boundary nuEffBf =
        turbModel.nuEff()().boundaryField();

    const volScalarField::Boundary nuBf =
        turbModel.nu()().boundaryField();

    const fvPatchList& patches = mesh_.boundary();

    volVectorField::Boundary& wallGradUMeanBf = wallGradUMean.boundaryFieldRef();

    IOobject UMeanheader
                (
                 "UMean",
                 mesh_.time().timeName(),
                 mesh_,
                 IOobject::MUST_READ
                );
    volVectorField UMean(UMeanheader, mesh_);

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
		wallGradUMeanBf[patchi] = -UMean.boundaryField()[patchi].snGrad();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallGradUMean::wallGradUMean
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log)
{
    volVectorField* wallGradUMeanPtr
    (
        new volVectorField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("0", dimless, Zero)
        )
    );

    mesh_.objectRegistry::store(wallGradUMeanPtr);

    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallGradUMean::~wallGradUMean()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallGradUMean::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    return true;
}


bool Foam::functionObjects::wallGradUMean::execute()
{
    volVectorField& wallGradUMean =
        mesh_.lookupObjectRef<volVectorField>(type());

    if (mesh_.foundObject<turbulenceModel>(turbulenceModel::propertiesName))
    {
        const turbulenceModel& model = mesh_.lookupObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

        calcWallGradUMean(model, wallGradUMean);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::wallGradUMean::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volVectorField& wallGradUMean =
        mesh_.lookupObject<volVectorField>(type());

    const volVectorField::Boundary& wallGradUMeanBf = wallGradUMean.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const vectorField& wallGradUMeanp = wallGradUMeanBf[patchi];

            const vector minWallGradUMean = gMin(wallGradUMeanp);
            const vector maxWallGradUMean = gMax(wallGradUMeanp);
            const vector avgWallGradUMean = gAverage(wallGradUMeanp);

            if (Pstream::master())
            {
                Log << "    patch " << patch.name()
                    << " y+ : min = " << minWallGradUMean << ", max = " << maxWallGradUMean
                    << ", average = " << avgWallGradUMean << nl;

                writeTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minWallGradUMean
                    << token::TAB << maxWallGradUMean
                    << token::TAB << avgWallGradUMean
                    << endl;
            }
        }
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
