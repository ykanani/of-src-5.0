/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "turbulentInflowProfileFvPatchField.H"
#include "clock.H"
#include "vectorList.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
turbulentInflowProfileFvPatchField<Type>::turbulentInflowProfileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    	fixedValueFvPatchField<Type>(p, iF),
        pointsDict_
                        (
                            IOobject
                            (
                                "points",
                                this->db().time().constant(),//timeName(this->db().time().timeOutputValue()-this->db().time().deltaTValue()),
                                this->db(),
                                IOobject::MUST_READ,
                                IOobject::AUTO_WRITE
                            )
                        ),

	ranGen_(clock().getTime()),
        n_(100),                        		//number of eddies
        sigmaMin_(0.5,0.5,0.5),         		//min length scale in all direction
        sigmaMax_(0.5,0.5,0.5),     	    	//max length scale in all direction     
        intensityMin_(1,1,1),  		//min intensity in all direction                
        intensityMax_(1,1,1),	   	//max intensity in all direction                
        Uinf_(1,0,0),   			//mean velocity
        L_(1),          				//eddy box length in flow direction
	referenceField_(p.size()),			//reference filed
        curTimeIndex_(-1)       			// time index  

{}


template<class Type>
turbulentInflowProfileFvPatchField<Type>::turbulentInflowProfileFvPatchField
(
    const turbulentInflowProfileFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
	fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
        pointsDict_(ptf.pointsDict_),
	ranGen_(clock().getTime()),
        n_(ptf.n_),
        sigmaMin_(ptf.sigmaMin_),
        sigmaMax_(ptf.sigmaMax_),
        intensityMin_(ptf.intensityMin_),
        intensityMax_(ptf.intensityMax_),
        Uinf_(ptf.Uinf_),
        L_(ptf.L_),
	referenceField_(ptf.referenceField_, mapper),
        curTimeIndex_(-1)
{}


template<class Type>
turbulentInflowProfileFvPatchField<Type>::turbulentInflowProfileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    	fixedValueFvPatchField<Type>(p, iF),
        pointsDict_
                        (
                            IOobject
                            (
                                "points",
                                this->db().time().constant(),//timeName(this->db().time().timeOutputValue()-this->db().time().deltaTValue()),
                                this->db(),
                                IOobject::MUST_READ,
                                IOobject::AUTO_WRITE
                            )
                        ),

	ranGen_(clock().getTime()),
        n_(readScalar(dict.lookup("n"))),
        sigmaMin_(dict.lookup("sigmaMin")),
        sigmaMax_(dict.lookup("sigmaMax")),
        intensityMin_(dict.lookup("intensityMin")),
        intensityMax_(dict.lookup("intensityMax")),
        Uinf_(dict.lookup("Uinf")),
        L_(readScalar(dict.lookup("L"))),
	referenceField_("referenceField", dict, p.size()),
        curTimeIndex_(-1)


{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==(referenceField_);
    }
}


template<class Type>
turbulentInflowProfileFvPatchField<Type>::turbulentInflowProfileFvPatchField
(
    const turbulentInflowProfileFvPatchField<Type>& ptf
)
:
    	fixedValueFvPatchField<Type>(ptf),
	pointsDict_(ptf.pointsDict_),
	ranGen_(ptf.ranGen_),
        n_(ptf.n_),
        sigmaMin_(ptf.sigmaMin_),
        sigmaMax_(ptf.sigmaMax_),
        intensityMin_(ptf.intensityMin_),
        intensityMax_(ptf.intensityMax_),
        Uinf_(ptf.Uinf_),
        L_(ptf.L_),
	referenceField_(ptf.referenceField_),
        curTimeIndex_(-1)
{}


template<class Type>
turbulentInflowProfileFvPatchField<Type>::turbulentInflowProfileFvPatchField
(
    const turbulentInflowProfileFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
	pointsDict_(ptf.pointsDict_),
	ranGen_(ptf.ranGen_),
        n_(ptf.n_),
        sigmaMin_(ptf.sigmaMin_),
        sigmaMax_(ptf.sigmaMax_),
        intensityMin_(ptf.intensityMin_),
        intensityMax_(ptf.intensityMax_),
        Uinf_(ptf.Uinf_),
        L_(ptf.L_),
        referenceField_(ptf.referenceField_),
        curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void turbulentInflowProfileFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    referenceField_.autoMap(m);
}


template<class Type>
void turbulentInflowProfileFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const turbulentInflowProfileFvPatchField<Type>& tiptf =
        refCast<const turbulentInflowProfileFvPatchField<Type> >(ptf);

    referenceField_.rmap(tiptf.referenceField_, addr);
}


template<class Type>
void turbulentInflowProfileFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
	
	Field<Type>& patchField = *this;
	Info << "in update coef" << endl;
	//reading and updating eddy positions and calculating patch velocity if it needs an update
	if (curTimeIndex_ != this->db().time().timeIndex()      )
        {

	// defining eddy data file
	
	// defining eddy data file 
                // getting bound box for patch
                boundBox bb(this->patch().patch().localPoints(), true);
                vector startPosition(bb.min()[0]-L_,bb.min()[1],bb.min()[2]);
//                vector endPosition(bb.min()[0],bb.max()[1],bb.max()[2]);
                vector endPosition(bb.min()[0],0.0043,bb.max()[2]);

                Info << "start = " << startPosition << endl;
                Info << "end = " << endPosition << endl;
                // constructing requierd variables
                vectorField pp(n_);
                vectorField sigma(n_);
                vectorField intensity(n_);
                vectorField rndsign(n_);
                vectorField signedintensity(n_);
		Field<Type> signedintensityType(n_);
		vector  tempvector(0.5,0.5,0.5);
                vector  tempsign(1,1,1);
                vector  unit(1,1,1);
		vector  tempi(0.01,0.01,0.01);
                //vector  tempvector(pTraits<Type>::one*0.5);
                //vector  tempsign(pTraits<Type>::one);
                //vector  unit(pTraits<Type>::one);
		//char&	typeName(pTraits<Type>::typeName);
		//genrating random positions for the first time
                bool isFirst(readBool(pointsDict_.lookup("isFirst")));
                if (isFirst)
                {

                        Info << "first time!!!" << endl;
                        //generating eddy positions
                        //scalar nEddy(readScalar(pointsDict_.lookup("nEddy")));


                        vector a(this->db().time().timeOutputValue(),1,1);      //dummy
                        forAll ( pp,i )
                        {
                                pp[i] = ranGen_.position(startPosition,endPosition);
                                //ranGen_.randomise(sigma[i]);
				//sigma[i] = sigma[i]*(sigmaMax_-sigmaMin_)+pTraits<Type>::one*sigmaMin_;
                                //ranGen_.randomise(intensity[i]);
				//intensity[i] = intensity[i]*(intensityMax_-intensityMin_)+pTraits<Type>::one*intensityMin_;
                                //ranGen_.randomise(tempsign);
                                //tempsign = tempsign - tempvector;

				//rndsign[i]=sign(tempsign);
                                sigma[i] = ranGen_.position(sigmaMin_,sigmaMax_);
                                intensity[i] = (intensityMax_-tempi)*(1-pp[i][1]/0.0043)+tempi;
                                ranGen_.randomise(tempsign);
                                tempsign = tempsign - tempvector;
                                rndsign[i][0]=sign(tempsign[0]);
                                rndsign[i][1]=sign(tempsign[1]);
                                rndsign[i][2]=sign(tempsign[2]);


                        }

                        pointsDict_.set("points",pp);
                        pointsDict_.set("sigma",sigma);
                        pointsDict_.set("intensity",intensity);
                        pointsDict_.set("rndsign",rndsign);
                        pointsDict_.set("time",this->db().time().timeOutputValue());
                        pointsDict_.set("isFirst",false);
                        pointsDict_.Foam::regIOobject::write();

                        //Info << "positions = " << pp << endl;
                        //Info << "sigma = " << sigma << endl;
                        //Info << "intensity = " << intensity << endl;
                        //Info << "rndsign = " << rndsign << endl;

                }
		else
                {
                        Info << "NOT first time!!!" << endl;
                        //reading positions             
                        pp              = vectorList(pointsDict_.lookup("points"));
                        sigma           =vectorList(pointsDict_.lookup("sigma"));
                        intensity       =vectorList(pointsDict_.lookup("intensity"));
                        rndsign         =vectorList(pointsDict_.lookup("rndsign"));

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
					pp[i][0]=startPosition[0];
                                        sigma[i] = ranGen_.position(sigmaMin_,sigmaMax_);
                                	intensity[i] = (intensityMax_-tempi)*(1-pp[i][1]/0.0043)+tempi;
                              		ranGen_.randomise(tempsign);
                                	tempsign = tempsign - tempvector;
                                	rndsign[i][0]=sign(tempsign[0]);
                                	rndsign[i][1]=sign(tempsign[1]);
                                	rndsign[i][2]=sign(tempsign[2]);


		
                                }


                        }



                        pointsDict_.set("points",pp);
                        pointsDict_.set("sigma",sigma);
                        pointsDict_.set("intensity",intensity);
                        pointsDict_.set("rndsign",rndsign);

                        pointsDict_.set("time",this->db().time().timeOutputValue());
                        pointsDict_.Foam::regIOobject::write();
                        //Info << "points = " << pp << endl;
                }



                //Info << "out of else " << endl;
                //Info << "positions = " << pp << endl;
                //Info << "sigma = " << sigma << endl;
                //Info << "intensity = " << intensity << endl;
                //Info << "rndsign = " << rndsign << endl;


                scalar Vb((max(pp + sigma)[0]-min(pp - sigma)[0])*(endPosition[1]-bb.min()[1])*(bb.max()[2]-bb.min()[2]));
                //Info << "maximum x = " << max(pp + sigma)[0] << endl;         
                //Info << "minimum x = " << min(pp - sigma)[0] << endl; 
                //Info << "Vb      = " << Vb << endl;
                tensor a(1,0,0,0,1,0,0,0,1);

                //calculating patch velocity based on eddies current poistion
                //vectorField& patchField = *this;
                const vectorField& c = this->patch().Cf();
                vector dx(0,0,0);
                scalar f(0);
		Type test;
	        vector m(0,0,0);

	        scalarField Coeff = this->patch().Cf() & m;

                forAll ( c,facei )
                {
                        patchField[facei] = pTraits<Type>::zero;//vector(0,0,0);
                        forAll(pp,i)
                        {
                                dx = c[facei]-pp[i];
                                dx[0]=dx[0]/sigma[i][0];
                                dx[1]=dx[1]/sigma[i][1];
                                dx[2]=dx[2]/sigma[i][2];


                                f=1;
                                for(int j=0; j<3; j=j+1)
                                {
                                        if (mag(dx[j]) >= 1 )
                                        {
                                                f = 0;
                                        }
                                        else
                                        {
                                                //fout = fout * cos(pi/2*x(i))^2;
                                                //f = f * sqrt(scalar(5))*exp(-5*mag(dx[j]));  //1
                                                f = f * sqrt(scalar(1.5))*(1-mag(dx[j]));      //T21
                                                //f = f * sqrt(scalar(3.0)/sqrt(scalar(3.141592)))*exp(-9/2*sqr(dx[j])); //T22
                                                //fout = fout * 2*exp(-9*x(i)^2/2);
                                        }
                                        signedintensity[i][j] = rndsign[i][j] * intensity[i][j];
                                }
                                //patchField[facei] = unit;
				
				
				signedintensityType[i] = simplify(signedintensity[i], pTraits<Type>::one);
				patchField[facei] = patchField[facei] +
                                                                        (
                                                                        (signedintensityType[i])*mag(Uinf_)
                                                                        /sqrt(sigma[i][0]*sigma[i][1]*sigma[i][2])
                                                                        * f
                                                                        * sqrt(Vb) / sqrt(n_)
                                                                        );

                                //scalar bb(0);
                                //bb = vector(sigma[i][0],0,0) & (vector(0,sigma[i][1],0) ^ vector(0,0,sigma[i][2]));//sigma[i] & vector(1,1,1);//a & signedintensity[i];
                                //Info << "result = " << sigma[i] << endl; //   /(sigma(1)*sigma(2)*sigma(3))* sqrt(Vb) *f((evalPosition - eddyPosition)./sigma);       

                        }
                        //adding mean velocity
                        scalar yOverDelta = ( c[facei][1]  )/0.0043 ;

                                    if (yOverDelta > 1.0)
                                            {
                                               Coeff[facei]= 1;

                                            }
                                    else
                                            {

                                                Coeff[facei]= 0.757* pow((c[facei][1]/0.0007), (1.0/6.7));

                                                }
 
			patchField[facei] = patchField[facei] + simplify(Coeff[facei]*Uinf_,pTraits<Type>::one);


                        //Info << "c[facei] =  " << c[facei] << endl;
                        //Info << "facei    =  " << facei << endl;
                        //Info << "result2  =  " << patchField[facei] << endl;
                }

                curTimeIndex_ = this->db().time().timeIndex();
        }
        /*
        else
        {
                        Info << "Do not change positions,calculate velocity based on existing positions " << endl;
                        //reading positions             
                        pp              =vectorList(pointsDict.lookup("points"));
                        sigma           =vectorList(pointsDict.lookup("sigma"));
                        intensity       =vectorList(pointsDict.lookup("intensity"));
                        rndsign         =vectorList(pointsDict.lookup("rndsign"));
                        //Info << "in else " << endl;
                        //Info << "positions = " << pp << endl;
                        //Info << "sigma = " << sigma << endl;
                        //Info << "intensity = " << intensity << endl;
                        //Info << "rndsign = " << rndsign << endl;
        }
        */















fixedValueFvPatchField<Type>::updateCoeffs();
}



template<class Type>
void turbulentInflowProfileFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("n") << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaMin") << sigmaMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaMax") << sigmaMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("intensityMin") << intensityMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("intensityMax") << intensityMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("Uinf") << Uinf_ << token::END_STATEMENT << nl;
    os.writeKeyword("L") << L_ << token::END_STATEMENT << nl;
    referenceField_.writeEntry("referenceField", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
