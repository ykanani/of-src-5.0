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
    Foam::functionObjects::wallGradUMean

Group
    grpFieldFunctionObjects

Description
    Evaluates and outputs turbulence y+ for models. Values written to
    time directories as field 'wallGradUMean'.

    Example of function object specification:
    \verbatim
    wallGradUMean1
    {
        type        wallGradUMean;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: wallGradUMean           | yes        |
    \endtable

Note
    Writing field 'wallGradUMean' is done by default, but it can be overridden by
    defining an empty \c objects list. For details see writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    wallGradUMean.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wallGradUMean_H
#define functionObjects_wallGradUMean_H

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
                          Class wallGradUMean Declaration
\*---------------------------------------------------------------------------*/

class wallGradUMean
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{
    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate y+
        void calcWallGradUMean
        (
            const turbulenceModel& turbModel,
            volVectorField& wallGradUMean
        );

        //- Disallow default bitwise copy construct
        wallGradUMean(const wallGradUMean&);

        //- Disallow default bitwise assignment
        void operator=(const wallGradUMean&);


public:

    //- Runtime type information
    TypeName("wallGradUMean");


    // Constructors

        //- Construct from Time and dictionary
        wallGradUMean
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~wallGradUMean();


    // Member Functions

        //- Read the wallGradUMean data
        virtual bool read(const dictionary&);

        //- Calculate the wallGradUMean field
        virtual bool execute();

        //- Write the wallGradUMean field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
