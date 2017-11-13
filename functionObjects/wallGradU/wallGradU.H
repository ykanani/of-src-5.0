/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::functionObjects::wallGradU

Group
    grpFieldFunctionObjects

Description
    Evaluates and outputs turbulence y+ for models. Values written to
    time directories as field 'wallGradU'.

    Example of function object specification:
    \verbatim
    wallGradU1
    {
        type        wallGradU;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: wallGradU           | yes        |
    \endtable

Note
    Writing field 'wallGradU' is done by default, but it can be overridden by
    defining an empty \c objects list. For details see writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    wallGradU.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wallGradU_H
#define functionObjects_wallGradU_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class turbulenceModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class wallGradU Declaration
\*---------------------------------------------------------------------------*/

class wallGradU
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{
    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate y+
        void calcWallGradU
        (
            const turbulenceModel& turbModel,
            volVectorField& wallGradU
        );

        //- Disallow default bitwise copy construct
        wallGradU(const wallGradU&);

        //- Disallow default bitwise assignment
        void operator=(const wallGradU&);


public:

    //- Runtime type information
    TypeName("wallGradU");


    // Constructors

        //- Construct from Time and dictionary
        wallGradU
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~wallGradU();


    // Member Functions

        //- Read the wallGradU data
        virtual bool read(const dictionary&);

        //- Calculate the wallGradU field
        virtual bool execute();

        //- Write the wallGradU field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
