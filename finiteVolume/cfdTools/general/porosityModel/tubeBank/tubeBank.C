/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "tubeBank.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(tubeBank, 0);
        addToRunTimeSelectionTable(porosityModel, tubeBank, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::tubeBank::tubeBank
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),
    sl_(readScalar(coeffs_.lookup("sl"))),
    lp_(readScalar(coeffs_.lookup("lp"))),
    //mC1_(readScalar(coeffs_.lookup("mC1"))),
    st_(readScalar(coeffs_.lookup("st"))),
    Dh_(readScalar(coeffs_.lookup("Dh"))),
    d_(readScalar(coeffs_.lookup("d"))),
    //C1_(readScalar(coeffs_.lookup("C1"))),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho"))

{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::tubeBank::~tubeBank()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::tubeBank::calcTransformModelData()
{
    // nothing to be transformed
}


void Foam::porosityModels::tubeBank::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    const scalarField& V = mesh_.V();
    Info << "calculating force..." << endl;
    apply(Udiag, V, rho, U);

    force = Udiag*U;
}


void Foam::porosityModels::tubeBank::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        apply(Udiag, V, rho, U);
        
    }
    else
    {
        apply(Udiag, V, geometricOneField(), U);
    }
}


void Foam::porosityModels::tubeBank::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
 
    apply(Udiag, V, rho, U);
}


void Foam::porosityModels::tubeBank::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        apply(AU, rho, U);
    }
    else
    {
        apply(AU, geometricOneField(), U);
    }
}


bool Foam::porosityModels::tubeBank::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}
// ************************************************************************* //
