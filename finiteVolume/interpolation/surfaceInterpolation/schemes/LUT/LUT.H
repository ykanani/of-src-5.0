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
    Foam::LUT

Description
    LUT: Linear-upwind transport.

    Interpolation scheme class derived from linearUpwind which returns blended
    linear/linear-upwind weighting factors and also applies a explicit
    gradient-based correction obtained from the linearUpwind scheme.  The
    blending-factor is set to 0.75 linear which optimises the balance between
    accuracy and stability on a range of LES cases with a range of mesh quality.

SourceFiles
    LUT.C

\*---------------------------------------------------------------------------*/

#ifndef LUT_H
#define LUT_H

#include "linearUpwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LUT Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LUT
:
    public linearUpwind<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        LUT(const LUT&);

	//- Percentage of linear scheme
	scalar cl_;
        //- Disallow default bitwise assignment
        void operator=(const LUT&);


public:

    //- Runtime type information
    TypeName("LUT");


    // Constructors

        //- Construct from mesh and Istream
        LUT
        (
            const fvMesh& mesh,
            Istream& schemeData
        )
        :
            linearUpwind<Type>(mesh, schemeData),
	    cl_(1.0)
        {}

        //- Construct from mesh, faceFlux and Istream
        LUT
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            linearUpwind<Type>(mesh, faceFlux, schemeData),
	    cl_(readScalar(schemeData))
        {
		//WarningInFunction << "cl_:" << cl_ << endl;
	}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return
                cl_*this->mesh().surfaceInterpolation::weights()
              + (1-cl_)*linearUpwind<Type>::weights();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return (1-cl_)*linearUpwind<Type>::correction(vf);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
