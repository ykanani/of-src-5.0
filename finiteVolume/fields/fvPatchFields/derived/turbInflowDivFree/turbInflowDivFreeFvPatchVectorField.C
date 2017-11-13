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

#include "turbInflowDivFreeFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "vectorList.H"
#include "clock.H"
#include "math.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbInflowDivFreeFvPatchVectorField::turbInflowDivFreeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
	fixedValueFvPatchVectorField(p, iF),
	ranGen_(clock().getTime()),
	n_(100),			//number of eddies
	sigmaMin_(0.1,0.1,0.1),		//min length scale in all direction
	sigmaMax_(0.5,0.5,0.5),		//max length scale in all direction	
	intensityMin_(1,1,1),	//min intensity in all direction		
	intensityMax_(1,1,1),	//max intensity in all direction		
	Uinf_(1,0,0),	//mean velocity
	L_(0),		//eddy box length in flow direction
	curTimeIndex_(-1)	// time index	
    
{}


turbInflowDivFreeFvPatchVectorField::turbInflowDivFreeFvPatchVectorField
(
    const turbInflowDivFreeFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    	fixedValueFvPatchVectorField(ptf, p, iF, mapper),
	ranGen_(ptf.ranGen_),
   	n_(ptf.n_),			
	sigmaMin_(ptf.sigmaMin_),		
	sigmaMax_(ptf.sigmaMax_),		
	intensityMin_(ptf.intensityMin_),		
	intensityMax_(ptf.intensityMax_),
	Uinf_(ptf.Uinf_),			
	L_(ptf.L_),
	curTimeIndex_(-1)
{}


turbInflowDivFreeFvPatchVectorField::turbInflowDivFreeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
	fixedValueFvPatchVectorField(p, iF),
	ranGen_(clock().getTime()),
	n_(readScalar(dict.lookup("n"))),			
	sigmaMin_(dict.lookup("sigmaMin")),		
	sigmaMax_(dict.lookup("sigmaMax")),		
	intensityMin_(dict.lookup("intensityMin")),		
	intensityMax_(dict.lookup("intensityMax")),
	Uinf_(dict.lookup("Uinf")),
	L_(readScalar(dict.lookup("L"))),
	curTimeIndex_(-1)

{}


turbInflowDivFreeFvPatchVectorField::turbInflowDivFreeFvPatchVectorField
(
    const turbInflowDivFreeFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
	fixedValueFvPatchVectorField(fcvpvf, iF),
	ranGen_(fcvpvf.ranGen_),
	n_(fcvpvf.n_),			
	sigmaMin_(fcvpvf.sigmaMin_),		
	sigmaMax_(fcvpvf.sigmaMax_),		
	intensityMin_(fcvpvf.intensityMin_),		
	intensityMax_(fcvpvf.intensityMax_),
	Uinf_(fcvpvf.Uinf_),			
	L_(fcvpvf.L_),
	curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbInflowDivFreeFvPatchVectorField::updateCoeffs()
{
	if (updated())
	{
		return;
	}


	////////////new

	//reading and updating eddy positions and calculating patch velocity if it needs an update
	if (curTimeIndex_ != this->db().time().timeIndex()	)
	{
		
		// defining eddy data file 
		IOdictionary  pointsDict
			(
			    IOobject
			    (
				"points",
				this->db().time().constant(),//timeName(this->db().time().timeOutputValue()-this->db().time().deltaTValue()),
				this->db(),
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE,
				false
			    )
			);

	
		// getting bound box for patch
		boundBox bb(patch().patch().localPoints(), true);
		vector startPosition(bb.min()[0]-L_,bb.min()[1],bb.min()[2]);	 
		vector endPosition(bb.min()[0],bb.max()[1],bb.max()[2]);

		Info << "start = " << startPosition << endl;
		Info << "end = " << endPosition << endl;
		// constructing requierd variables
		vectorField pp(n_);
		vectorField sigma(n_);
		vectorField intensity(n_);
		vectorField rndsign(n_);
		vectorField signedintensity(n_);
		vector	tempvector(0.5,0.5,0.5);
		vector  tempsign(1,1,1);
		vector  unit(1,1,1);





		//genrating random positions for the first time
		bool isFirst(readBool(pointsDict.lookup("isFirst")));
		if (isFirst)
		{	
		
			Info << "first time!!!" << endl;		
			//generating eddy positions
			//scalar nEddy(readScalar(pointsDict.lookup("nEddy")));
			

			vector a(this->db().time().timeOutputValue(),1,1);	//dummy
			forAll ( pp,i )
			{
				pp[i] = ranGen_.position(startPosition,endPosition);
				sigma[i] = ranGen_.position(sigmaMin_,sigmaMax_);
				intensity[i] = ranGen_.position(intensityMin_,intensityMax_);
				ranGen_.randomise(tempsign);
				tempsign = tempsign - tempvector;
				rndsign[i][0]=sign(tempsign[0]);
				rndsign[i][1]=sign(tempsign[1]);
				rndsign[i][2]=sign(tempsign[2]);			

			}
		
			pointsDict.set("points",pp);
			pointsDict.set("sigma",sigma);
			pointsDict.set("intensity",intensity);
			pointsDict.set("rndsign",rndsign);
			pointsDict.set("time",this->db().time().timeOutputValue());
			pointsDict.set("isFirst",false);
			pointsDict.Foam::regIOobject::write();
		
			//Info << "positions = " << pp << endl;
			//Info << "sigma = " << sigma << endl;
			//Info << "intensity = " << intensity << endl;
			//Info << "rndsign = " << rndsign << endl;

		}
		else
		{
			Info << "NOT first time!!!" << endl;		
			//reading positions		
			pp		= vectorList(pointsDict.lookup("points"));
			sigma		=vectorList(pointsDict.lookup("sigma"));
			intensity	=vectorList(pointsDict.lookup("intensity"));
			rndsign 	=vectorList(pointsDict.lookup("rndsign"));

			//convecting eddies with mean velocity
			vector a(this->db().time().timeOutputValue(),1,1);   //dummy
			forAll ( pp,i )
			{
				pp[i]=pp[i]+Uinf_*this->db().time().deltaTValue();
				//checking if eddies convected outside of the box
				if ( pp[i][0] - sigma[i][0] > bb.min()[0] )
				{
					//generating new eddy
					//Info << "generating new eddy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1" << endl;
					pp[i] = ranGen_.position(startPosition,endPosition);
					sigma[i] = ranGen_.position(sigmaMin_,sigmaMax_);
					intensity[i] = ranGen_.position(intensityMin_,intensityMax_);
					ranGen_.randomise(tempsign);
					tempsign = tempsign - tempvector;
					rndsign[i][0]=sign(tempsign[0]);
					rndsign[i][1]=sign(tempsign[1]);
					rndsign[i][2]=sign(tempsign[2]);
				}
					 				

			}

			
			
                	pointsDict.set("points",pp);
			pointsDict.set("sigma",sigma);
			pointsDict.set("intensity",intensity);
			pointsDict.set("rndsign",rndsign);

			pointsDict.set("time",this->db().time().timeOutputValue());
			pointsDict.Foam::regIOobject::write();
			//Info << "points = " << pp << endl;
		}
		
		//Info << "out of else " << endl;
		//Info << "positions = " << pp << endl;
		//Info << "sigma = " << sigma << endl;
		//Info << "intensity = " << intensity << endl;
		//Info << "rndsign = " << rndsign << endl;
		
	
		scalar Vb((max(pp + sigma)[0]-min(pp - sigma)[0])*(bb.max()[1]-bb.min()[1])*(bb.max()[2]-bb.min()[2]));
		//Info << "maximum x = " << max(pp + sigma)[0] << endl;		
		//Info << "minimum x = " << min(pp - sigma)[0] << endl;	
		//Info << "Vb	   = " << Vb << endl;
		tensor a(1,0,0,0,1,0,0,0,1);
		
		//calculating patch velocity based on eddies current poistion
		vectorField& patchField = *this;
		const vectorField& c = patch().Cf();
		vector dx(0,0,0);
		scalar q(0);
		scalar rsigma(0);
		forAll ( c,facei ) 
		{
			patchField[facei] = vector(0,0,0);
			forAll(pp,i)
			{
				//calculating distance to eddy (r)
				dx = c[facei]-pp[i];
				dx[0]=dx[0]/sigma[i][0];
				dx[1]=dx[1]/sigma[i][1];
				dx[2]=dx[2]/sigma[i][2];
				
				rsigma = mag(dx)/mag(sigma[i]);
				Info << "rsigma" << rsigma << endl;

				q=1;
				if ( rsigma > 1 )
				{
					q = 0;
				}
				else
				{
					//fout = fout * cos(pi/2*x(i))^2;
						q =  	
						sqrt(intensity[i][0]*mag(Uinf_))
						*sqrt(16*Vb/( 15*3.14*pow(mag(sigma[i]),3)))*sqr(sin(3.14*rsigma))
						/sqr(rsigma)
						/mag(sigma[i]) ;
					//fout = fout * sqrt(3/2)*(1-abs(x(i)^.1));
					//fout = fout * 2*exp(-9*x(i)^2/2);
				}
					//signedintensity[i][j] = rndsign[i][j] * intensity[i][j];					
				//patchField[facei] = unit;
					patchField[facei] = patchField[facei] +
									(
									q * (dx ^ rndsign[i])
									);
				
				//scalar bb(0);
				//bb = vector(sigma[i][0],0,0) & (vector(0,sigma[i][1],0) ^ vector(0,0,sigma[i][2]));//sigma[i] & vector(1,1,1);//a & signedintensity[i];
				//Info << "result = " << sigma[i] << endl; //   /(sigma(1)*sigma(2)*sigma(3))* sqrt(Vb) *f((evalPosition - eddyPosition)./sigma);	
				
			}
			//adding mean velocity
						patchField[facei] = patchField[facei] * sqrt(1/n_) + Uinf_;
			//Info << "c[facei] =  " << c[facei] << endl;
			//Info << "facei    =  " << facei << endl;
			//Info << "result2  =  " << patchField[facei] << endl;
		}
		Info << "patch average" <<  average(patchField) << endl;

		curTimeIndex_ = this->db().time().timeIndex();
	}
	/*
	else
	{
			Info << "Do not change positions,calculate velocity based on existing positions " << endl;
			//reading positions		
			pp		=vectorList(pointsDict.lookup("points"));
			sigma		=vectorList(pointsDict.lookup("sigma"));
			intensity	=vectorList(pointsDict.lookup("intensity"));
			rndsign 	=vectorList(pointsDict.lookup("rndsign"));
			//Info << "in else " << endl;
			//Info << "positions = " << pp << endl;
			//Info << "sigma = " << sigma << endl;
			//Info << "intensity = " << intensity << endl;
			//Info << "rndsign = " << rndsign << endl;
	}
	*/
	



	fixedValueFvPatchVectorField::updateCoeffs();









/////end new

/*
//Info << "delta_ = " << delta_ << endl;
	//Info << "Uzero_ = " << Uzero_ << endl;
	//Info << "power_ = " << power_ << endl;
    // initialise a zero−filled field of the same 
    // size as the list of face centres of this patch 
	//vector n(1,0,0);
	vector m(0,0,0);
	vector n(1,0,0);
	scalarField Coeff = patch().Cf() & m;

	//boundBox bb(patch().patch().localPoints(), true);
	scalar y0 = bb.min()[1];
	//Info << "bb.min = " << bb.min() << endl;
	//Info << "y0 = " << y0 << endl;
    // go over each face 
    forAll ( c,facei ) 
    {
	// define the non−dimensional distancce to 
        // establish where in the BL this face sits
        //Info << "c" << c[facei][1] << endl;
	//scalar yOverDelta = ( c[facei][1] - y0 )/delta_ ;
	
	//Info << "yoverdelta = " << yOverDelta << endl;

	// if we are outside of BL, assign free-streem speed
    if (yOverDelta > 1.0)
    {
       Coeff[facei]= mag(Uzero_);
	
    }
    else
    {
        
        Coeff[facei]= mag(Uzero_);// * pow(yOverDelta, power_);
	
    }
	
    }

	vectorField::operator=(n*Coeff);
*/
   
}


// Write
void turbInflowDivFreeFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaMin")
        << sigmaMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaMax")
        << sigmaMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("intensityMin")
        << intensityMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("intensityMax")
        << intensityMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("Uinf")
        << Uinf_ << token::END_STATEMENT << nl;
    os.writeKeyword("L")
        << L_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, turbInflowDivFreeFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
