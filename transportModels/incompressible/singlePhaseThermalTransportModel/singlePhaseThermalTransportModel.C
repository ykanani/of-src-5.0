/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "singlePhaseThermalTransportModel.H"
#include "viscosityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singlePhaseThermalTransportModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseThermalTransportModel::singlePhaseThermalTransportModel
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    viscosityModelPtr_(viscosityModel::New("nu", *this, U, phi))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singlePhaseThermalTransportModel::~singlePhaseThermalTransportModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::singlePhaseThermalTransportModel::nu() const
{
    return viscosityModelPtr_->nu();
}


Foam::tmp<Foam::scalarField>
Foam::singlePhaseThermalTransportModel::nu(const label patchi) const
{
    return viscosityModelPtr_->nu(patchi);
}

Foam::tmp<Foam::volScalarField>
Foam::singlePhaseThermalTransportModel::alpha() const
{
    return viscosityModelPtr_->alpha();
}


Foam::tmp<Foam::scalarField>
Foam::singlePhaseThermalTransportModel::alpha(const label patchi) const
{
    return viscosityModelPtr_->alpha(patchi);
}


void Foam::singlePhaseThermalTransportModel::correct()
{
    viscosityModelPtr_->correct();
}


bool Foam::singlePhaseThermalTransportModel::read()
{
    if (regIOobject::read())
    {
        return viscosityModelPtr_->read(*this);
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
