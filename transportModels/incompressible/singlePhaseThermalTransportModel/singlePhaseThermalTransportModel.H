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
    Foam::singlePhaseThermalTransportModel

Description
    A simple single-phase transport model based on viscosityModel.

    Used by the incompressible single-phase solvers like simpleFoam,
    turbFoam etc.

SourceFiles
    singlePhaseThermalTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseThermalTransportModel_H
#define singlePhaseThermalTransportModel_H

#include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class viscosityModel;

/*---------------------------------------------------------------------------*\
                Class singlePhaseThermalTransportModel Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseThermalTransportModel
:
    public IOdictionary,
    public transportModel
{
    // Private Data

        autoPtr<viscosityModel> viscosityModelPtr_;


    // Private Member Functions

        //- Disallow copy construct
        singlePhaseThermalTransportModel(const singlePhaseThermalTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const singlePhaseThermalTransportModel&);


public:

    //- Runtime type information
    TypeName("singlePhaseThermalTransportModel");


    // Constructors

        //- Construct from components
        singlePhaseThermalTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~singlePhaseThermalTransportModel();


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const;

        //- Return the laminar viscosity
        virtual tmp<volScalarField> alpha() const;


        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> alpha(const label patchi) const;


        //- Correct the laminar viscosity
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
