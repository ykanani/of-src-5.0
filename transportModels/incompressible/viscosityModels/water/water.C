/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "water.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(water, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        water,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::water::calcNu() const
{
    const volScalarField& T= U_.mesh().lookupObject<volScalarField>("T1");
    const volScalarField Tnon = T/dimensionedScalar("",dimTemperature,1); 
    return (3.75447887045759E-14*pow(Tnon,4) - 5.0229976347443E-11*pow(Tnon,3)  + 2.52828128903637E-08*pow(Tnon,2)  - 5.68131312353983E-06*Tnon + 0.000481893250547309)*dimensionedScalar("",dimensionSet(0,2,-1,0,0),1);
//    return 2.414e-5*pow(10,247.8/(Tnon-140))*dimensionedScalar("",dimensionSet(0,2,-1,0,0),1);

    
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::water::calcAlpha() const
{
    const volScalarField& T= U_.mesh().lookupObject<volScalarField>("T1");
    const volScalarField Tnon = T/dimensionedScalar("",dimTemperature,1); 
    return  (-1.56206080872566E-16*pow(Tnon,4)+ 1.95710425241534E-13*pow(Tnon,3) - 9.39246289301658E-11*pow(Tnon,2) + 2.08284501467685E-08*Tnon - 1.66766616986267E-06)*dimensionedScalar("",dimensionSet(0,2,-1,0,0),1);

    
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::water::water
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
//    waterCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
   // k_("k", dimViscosity, waterCoeffs_),
   // n_("n", dimless, waterCoeffs_),
   // nuMin_("nuMin", dimViscosity, waterCoeffs_),
   // nuMax_("nuMax", dimViscosity, waterCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    ),
    alpha_
    (
        IOobject
        (
            "alpha",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcAlpha()
    )

{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::water::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

 //   waterCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

 //   waterCoeffs_.lookup("k") >> k_;
 //   waterCoeffs_.lookup("n") >> n_;
 //   waterCoeffs_.lookup("nuMin") >> nuMin_;
 //   waterCoeffs_.lookup("nuMax") >> nuMax_;

    return true;
}


// ************************************************************************* //
