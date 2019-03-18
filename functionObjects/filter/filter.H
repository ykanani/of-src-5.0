/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::functionObjects::filter

Group
    grpFieldFunctionObjects

Description
    Calculates the filterient of a field.

    The operation is limited to scalar and vector volume or surface fields, and
    the output is a volume vector or tensor field.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    filter.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_filter_H
#define functionObjects_filter_H

#include "fieldExpression.H"
#include "LESfilter.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class filter Declaration
\*---------------------------------------------------------------------------*/

class filter
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the magnitude of the field and register the result
        template<class Type>
        bool calcFilter();

        //- Calculate the filterient field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("filter");


    // Constructors

        //- Construct from Time and dictionary
        filter
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~filter();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "filterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
