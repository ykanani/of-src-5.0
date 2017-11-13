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

\*---------------------------------------------------------------------------*/

#include "powerLawVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

powerLawVelocityFvPatchVectorField::powerLawVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Uzero_(10,0,0),
    delta_(1),
    power_(1)
{}


powerLawVelocityFvPatchVectorField::powerLawVelocityFvPatchVectorField
(
    const powerLawVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Uzero_(ptf.Uzero_),
    delta_(ptf.delta_),
    power_(ptf.power_)
{}


powerLawVelocityFvPatchVectorField::powerLawVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    Uzero_(dict.lookup("Uzero")),
    delta_(readScalar(dict.lookup("delta"))),
    power_(readScalar(dict.lookup("power")))

{}


powerLawVelocityFvPatchVectorField::powerLawVelocityFvPatchVectorField
(
    const powerLawVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    Uzero_(fcvpvf.Uzero_),
    delta_(fcvpvf.delta_),
    power_(fcvpvf.power_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void powerLawVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	//Info << "delta_ = " << delta_ << endl;
	//Info << "Uzero_ = " << Uzero_ << endl;
	//Info << "power_ = " << power_ << endl;
    // initialise a zero−filled field of the same 
    // size as the list of face centres of this patch 
	//vector n(1,0,0);
	vector m(0,0,0);

	const vectorField& c = patch().Cf();
	scalarField Coeff = patch().Cf() & m;

	boundBox bb(patch().patch().localPoints(), true);
	scalar y0 = bb.min()[1];
	//Info << "bb.min = " << bb.min() << endl;
	//Info << "y0 = " << y0 << endl;
    // go over each face 
    forAll ( c,facei ) 
    {
	// define the non−dimensional distancce to 
        // establish where in the BL this face sits
        //Info << "c" << c[facei][1] << endl;
	scalar yOverDelta = ( c[facei][1] - y0 )/delta_ ;
	
	//Info << "yoverdelta = " << yOverDelta << endl;

	// if we are outside of BL, assign free-streem speed
    if (yOverDelta > 1.0)
    {
       Coeff[facei]= mag(Uzero_);
	
    }
    else
    {
        
        Coeff[facei]= mag(Uzero_) * pow(yOverDelta, power_);
	
    }
	
    }

	vectorField::operator=(Uzero_/mag(Uzero_)*Coeff);
   
}


// Write
void powerLawVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Uzero")
        << Uzero_ << token::END_STATEMENT << nl;
    os.writeKeyword("delta")
        << delta_ << token::END_STATEMENT << nl;
    os.writeKeyword("power")
        << power_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, powerLawVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
