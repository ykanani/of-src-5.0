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

#include "wallGradU.H"
#include "turbulenceModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallGradU, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        wallGradU,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::wallGradU::writeFileHeader(const label i)
{
    writeHeader(file(), "y+ ()");

    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "average");
    file() << endl;
}


void Foam::functionObjects::wallGradU::calcWallGradU
(
    const turbulenceModel& turbModel,
    volVectorField& wallGradU
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

    volVectorField::Boundary& wallGradUBf = wallGradU.boundaryFieldRef();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            wallGradUBf[patchi] = -turbModel.U().boundaryField()[patchi].snGrad();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallGradU::wallGradU
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
    volVectorField* wallGradUPtr
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

    mesh_.objectRegistry::store(wallGradUPtr);

    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallGradU::~wallGradU()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallGradU::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    return true;
}


bool Foam::functionObjects::wallGradU::execute()
{
    volVectorField& wallGradU =
        mesh_.lookupObjectRef<volVectorField>(type());

    if (mesh_.foundObject<turbulenceModel>(turbulenceModel::propertiesName))
    {
        const turbulenceModel& model = mesh_.lookupObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

        calcWallGradU(model, wallGradU);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::wallGradU::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volVectorField& wallGradU =
        mesh_.lookupObject<volVectorField>(type());

    const volVectorField::Boundary& wallGradUBf = wallGradU.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const vectorField& wallGradUp = wallGradUBf[patchi];

            const vector minWallGradU = gMin(wallGradUp);
            const vector maxWallGradU = gMax(wallGradUp);
            const vector avgWallGradU = gAverage(wallGradUp);

            if (Pstream::master())
            {
                Log << "    patch " << patch.name()
                    << " y+ : min = " << minWallGradU << ", max = " << maxWallGradU
                    << ", average = " << avgWallGradU << nl;

                writeTime(file());
                file()
                    << token::TAB << patch.name()
                    << token::TAB << minWallGradU
                    << token::TAB << maxWallGradU
                    << token::TAB << avgWallGradU
                    << endl;
            }
        }
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
