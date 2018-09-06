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

Class
    Foam::viscosityModels::water

Description
     Standard power-law non-Newtonian viscosity model.

SourceFiles
    water.C

\*---------------------------------------------------------------------------*/

#ifndef water_H
#define water_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class water Declaration
\*---------------------------------------------------------------------------*/

class water
:
    public viscosityModel
{
    // Private data

//        dictionary waterCoeffs_;

//        dimensionedScalar k_;
//        dimensionedScalar n_;
//        dimensionedScalar nuMin_;
//        dimensionedScalar nuMax_;

        volScalarField nu_;
        volScalarField alpha_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;
        tmp<volScalarField> calcAlpha() const;


public:

    //- Runtime type information
    TypeName("water");


    // Constructors

        //- Construct from components
        water
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~water()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }
        //- Return the laminar thermal diffusivity
        virtual tmp<volScalarField> alpha() const
        {
            return alpha_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }
        //- Return the laminar thermal diffusivity for patch
        virtual tmp<scalarField> alpha(const label patchi) const
        {
            return alpha_.boundaryField()[patchi];
        }


        //- Correct the laminar viscosity
        virtual void correct()
        {
            nu_ = calcNu();
            alpha_ = calcAlpha();

        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
