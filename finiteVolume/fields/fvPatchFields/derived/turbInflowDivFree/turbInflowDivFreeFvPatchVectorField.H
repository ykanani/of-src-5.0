/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    parabolicVelocityFvPatchVectorField

Description
    Boundary condition specifies a parabolic velocity inlet profile
    (fixed value), given maximum velocity value (peak of the parabola),
    flow direction n and direction of the parabolic coordinate y

SourceFiles
    parabolicVelocityFvPatchVectorField.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef turbInflowDivFreeFvPatchVectorField_H
#define turbInflowDivFreeFvPatchVectorField_H

#include "fvPatchFields.H"
#include "Random.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class parabolicVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class turbInflowDivFreeFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
	
        Random	ranGen_;		//- Random number generator
	scalar	n_;		//number of eddies
	vector	sigmaMin_;		//min length scale in all direction
	vector	sigmaMax_;		//max length scale in all direction	

	vector	intensityMin_;	//min intensity in all direction		
	vector	intensityMax_;	//max intensity in all direction		


	vector  Uinf_;		//mean velocity	
        scalar	L_;		//eddy box length in flow direction
	label curTimeIndex_;	//- Current time index (used for updating)


public:

    //- Runtime type information
    TypeName("turbInflowDivFree");


    // Constructors

        //- Construct from patch and internal field
        turbInflowDivFreeFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbInflowDivFreeFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbInflowFvPatchVectorField
        //  onto a new patch
        turbInflowDivFreeFvPatchVectorField
        (
            const turbInflowDivFreeFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbInflowDivFreeFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbInflowDivFreeFvPatchVectorField
        (
            const turbInflowDivFreeFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new turbInflowDivFreeFvPatchVectorField(*this, iF)
            );
        }


    // Member functions
	
	//return mu,ber of eddies
	scalar&	n()
	{
		return n_;
	}

	//return min length scale in all direction		
	vector&	sigmaMin()
	{
		return sigmaMin_;
	}
	//return max length scale in all direction			
	vector&	sigmaMax()
	{
		return sigmaMax_;
	}
	//return min intensity in all direction		
	vector&	intensityMin()	
	{
		return intensityMin_;
	}
	//max intensity in all direction		
	vector&	intensityMax()
	{
		return intensityMax_;
	}
	//return eddy box length in flow direction
	scalar&	L()	
	{
		return L_;
	}
        //return mean velocity
	vector&	Uinf()	
	{
		return Uinf_;
	}
        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
