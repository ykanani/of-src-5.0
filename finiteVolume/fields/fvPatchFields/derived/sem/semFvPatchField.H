/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::semFvPatchField

Group
    grpInletBoundaryConditions

Description
    This boundary condition generates a fluctuating inlet condition by adding
    a random component to a reference (mean) field.

    \f[
        x_p = (1 - \alpha) x_p^{n-1} + \alpha (x_{ref} + s C_{RMS} x_{ref})
    \f]

    where

    \vartable
        x_p     | patch values
        x_{ref} | reference patch values
        n       | time level
        \alpha  | fraction of new random component added to previous time value
        C_{RMS} | RMS coefficient
        s       | fluctuation scale
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        fluctuationScale | RMS fluctuation scale (fraction of mean) | yes |
        referenceField | reference (mean) field | yes        |
        alpha | fraction of new random component added to previous| no| 0.1
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            sem;
        fluctuationScale 0.1;
        referenceField  uniform 10;
        alpha           0.1;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    semFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef semFvPatchField_H
#define semFvPatchField_H
#include "fvPatchFields.H"
#include "Random.H"
#include "pointToPointPlanarInterpolation.H"
#include "fixedValueFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class semFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class semFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data
	//- 2D interpolation (for 'planarInterpolation' mapMethod)
	autoPtr<pointToPointPlanarInterpolation> mapperPtr_;
	
	//- Points dictionary
	IOdictionary  pointsDict_;

	//- Input velocity and Reynolds stress tensor
	IOdictionary  statisticsDict_;

        //- Random number generator
        Random ranGen_;

	//- Number of eddies
	scalar	n_;

	//- Min Length scale
	vector	sigmaMin_;

	//- Max length scale
	vector	sigmaMax_;


        //- Mean velocity
        vector	Uinf_;
	
	//- friction coefficient Cf
	scalar	cf_;
	//- Eddy box length in flow direction
	scalar	L_;
	
        //- Reference field
        Field<Type> referenceField_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        vectorField pp_;
        vectorField sigma_;
        vectorField intensity_;
        vectorField rndsign_;
	scalar delta_;
	scalar fstu_;		//freestream turbulence
	scalar maxy_;


public:

    //- Runtime type information
    TypeName("sem");


    // Constructors

        //- Construct from patch and internal field
        semFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        semFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given semFvPatchField
        //  onto a new patch
        semFvPatchField
        (
            const semFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        semFvPatchField
        (
            const semFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new semFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        semFvPatchField
        (
            const semFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new semFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access
/*
            //- Return the fluctuation scale
            const Type& fluctuationScale() const
            {
                return fluctuationScale_;
            }

            //- Return reference to the fluctuation scale to allow adjustment
            Type& fluctuationScale()
            {
                return fluctuationScale_;
            }
*/
            //- Return the reference field
            const Field<Type>& referenceField() const
            {
                return referenceField_;
            }

            //- Return reference to the reference field to allow adjustment
            Field<Type>& referenceField()
            {
                return referenceField_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

	    vector sign2(vector D);
	    scalar sign2(scalar D)
		{
			return (D >= 0)? 1: -1;	
		}
            scalar simplify(vector I, scalar s)
                {
                        return I.x();
                }
            vector simplify(vector I, vector v)
                {
                        return I;
                }
            sphericalTensor simplify(vector I, sphericalTensor v)
                {
                        return v;
                }
            symmTensor  simplify(vector I, symmTensor v)
                {
                        return v;
                }
            tensor      simplify(vector I, tensor v)
                {
                        return v;
                }
     	    scalar dotp(tensor I, scalar s)
                {
                        return s;
                }
            vector dotp(tensor I, vector v)
                {
                        return I&v;
                }
            sphericalTensor dotp(tensor I, sphericalTensor v)
                {
                        return v;
                }
            symmTensor  dotp(tensor I, symmTensor v)
                {
                        return v;
                }
            tensor      dotp(tensor I, tensor v)
                {
                        return v;
                }



        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "semFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
