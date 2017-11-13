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

#include "turbulentInflowAmesFvPatchField.H"
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
turbulentInflowAmesFvPatchField<Type>::turbulentInflowAmesFvPatchField
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
        curTimeIndex_(-1),       			// time index 
	pp_(100),
        sigma_(100),
        intensity_(100),
        rndsign_(100),
	delta_(0.001),
	fstu_(0),
	maxy_(0.01)
 

{}


template<class Type>
turbulentInflowAmesFvPatchField<Type>::turbulentInflowAmesFvPatchField
(
    const turbulentInflowAmesFvPatchField<Type>& ptf,
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
        curTimeIndex_(-1),
        pp_(ptf.pp_),
        sigma_(ptf.sigma_),
        intensity_(ptf.intensity_),
        rndsign_(ptf.rndsign_),
	delta_(ptf.delta_),
	fstu_(ptf.fstu_),
	maxy_(ptf.maxy_)
{}


template<class Type>
turbulentInflowAmesFvPatchField<Type>::turbulentInflowAmesFvPatchField
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
        curTimeIndex_(-1),
        pp_(vectorList(pointsDict_.lookup("points"))),
        sigma_(vectorList(pointsDict_.lookup("sigma"))),
        intensity_(vectorList(pointsDict_.lookup("intensity"))),
        rndsign_(vectorList(pointsDict_.lookup("rndsign"))),
	delta_(readScalar(dict.lookup("delta"))),
	fstu_(readScalar(dict.lookup("fstu"))),
 	maxy_(readScalar(dict.lookup("maxy")))

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
turbulentInflowAmesFvPatchField<Type>::turbulentInflowAmesFvPatchField
(
    const turbulentInflowAmesFvPatchField<Type>& ptf
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
        curTimeIndex_(-1),
        pp_(ptf.pp_),
        sigma_(ptf.sigma_),
        intensity_(ptf.intensity_),
        rndsign_(ptf.rndsign_),
	delta_(ptf.delta_),
	fstu_(ptf.fstu_),
 	maxy_(ptf.maxy_)

{}


template<class Type>
turbulentInflowAmesFvPatchField<Type>::turbulentInflowAmesFvPatchField
(
    const turbulentInflowAmesFvPatchField<Type>& ptf,
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
        curTimeIndex_(-1),
        pp_(ptf.pp_),
        sigma_(ptf.sigma_),
        intensity_(ptf.intensity_),
        rndsign_(ptf.rndsign_),
	delta_(ptf.delta_),
	fstu_(ptf.fstu_),
	maxy_(ptf.maxy_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void turbulentInflowAmesFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    referenceField_.autoMap(m);
}


template<class Type>
void turbulentInflowAmesFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const turbulentInflowAmesFvPatchField<Type>& tiptf =
        refCast<const turbulentInflowAmesFvPatchField<Type> >(ptf);

    referenceField_.rmap(tiptf.referenceField_, addr);
}


template<class Type>
void turbulentInflowAmesFvPatchField<Type>::updateCoeffs()
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
		//experiments profile
		int n_inter;
		n_inter = 36;
  		scalar xp=0.0075005191;
		scalar yp=0.0546987521;
		scalar yf=0;
		scalar xf=0;
		scalar d = 0;	
                vector ppr(0,0,0);
		vector Uconv(0,0,0);
		// getting bound box for patch
                boundBox bb(this->patch().patch().localPoints(), true);
                vector bbmin(bb.min()[0]-L_,bb.min()[1],bb.min()[2]);
                vector bbmax(bb.min()[0],bb.max()[1],bb.max()[2]);

		vector p0(bb.max()[0],bb.min()[1],bb.min()[2]);
		vector p1(bb.min()[0],bb.max()[1],bb.max()[2]);
		vector p10(p1-p0);
		scalar angle = atan(float( (p0[0]-p1[0]) / (p1[1]-p0[1])  ));
		vector p10r(0,0,0);
		tensor rotate1(cos(-angle),-sin(-angle),0,sin(-angle),cos(-angle),0,0,0,1); 	// rotation from original to fake coordinate
		tensor rotate2(cos(angle),-sin(angle),0,sin(angle),cos(angle),0,0,0,1);		// rotation from fake to original coordinate

		p10r = rotate1 & p10;
		vector p1r(p0+p10r); 
		//Info << "p0 = " << p0 << endl;
		//Info << "p1 = " << p1 << endl;
		//Info << "angle = " << angle << endl;
                //Info << "p1r = " << p1r << endl;

                vector startPosition(p0[0]-L_,p0[1],p0[2]);
                vector endPosition(p0[0],p1r[1],p1r[2]);

                Info << "start = " << startPosition << endl;
                Info << "end = " << endPosition << endl;

                vectorField signedintensity(n_);
		Field<Type> signedintensityType(n_);
		vector  tempvector(0.5,0.5,0.5);
                vector  tempsign(1,1,1);
                vector  unit(1,1,1);
		vector  fstu(fstu_,fstu_,fstu_);

                bool isFirst(readBool(pointsDict_.lookup("isFirst")));
                if (isFirst)
                {

                        Info << "first time!!!" << endl;


                        vector a(this->db().time().timeOutputValue(),1,1);      //dummy
                        forAll ( pp_,i )
                        {
                                pp_[i] = ranGen_.position(startPosition,endPosition);	//generate point in fake coordinate
                               	//Info << "positions = " << pp_[i] << endl; 
                                sigma_[i] = ranGen_.position(sigmaMin_,sigmaMax_);
				intensity_[i] = fstu;
				if (pp_[i][1]<=delta_) //// no flactuation in the boundary layer
				{
                                	intensity_[i] = vector(0,0,0);
				}
				pp_[i] = (rotate2 & (pp_[i]-p0)) + p0;	//rotate it to the originat coordinate
                                ranGen_.randomise(tempsign);
                                tempsign = tempsign - tempvector;
                                rndsign_[i][0]=sign(tempsign[0]);
                                rndsign_[i][1]=sign(tempsign[1]);
                                rndsign_[i][2]=sign(tempsign[2]);


                        }

                        pointsDict_.set("points",pp_);
                        pointsDict_.set("sigma",sigma_);
                        pointsDict_.set("intensity",intensity_);
                        pointsDict_.set("rndsign",rndsign_);
                        pointsDict_.set("time",this->db().time().timeOutputValue());
                        pointsDict_.set("isFirst",false);
                        pointsDict_.Foam::regIOobject::write();

                        //Info << "positions = " << pp_ << endl;
                        //Info << "sigma = " << sigma_ << endl;
                        //Info << "intensity = " << intensity_ << endl;
                        //Info << "rndsign = " << rndsign_ << endl;

                }
		else
                {
                        Info << "NOT first time!!!" << endl;
                        //convecting eddies with mean velocity
                        vector a(this->db().time().timeOutputValue(),1,1);   //dummy
                        forAll ( pp_,i )
                        {
	                        xf = pp_[i][0]; // face i, x location
	                        yf = pp_[i][1]; // face i, y location
        	                d = sqrt( pow(xf-xp,2)+ pow(yf-yp,2) );

				if ( d < 0.002 )
				{
					Uconv[0]=1.785*d/0.002;
					Uconv[1]=5.635*d/0.002;
					Uconv[2]=0;

				}
				else
				{
					Uconv[0]=780.95*pow(d,3) - 431.07*pow(d,2) + 91.56*d + 1.6062;
					Uconv[1]=-175.24*pow(d,3) + 108.65*pow(d,2) - 32.449*d + 5.7001;
					Uconv[2]=0;

				}

                
                                pp_[i]=pp_[i]+Uconv*this->db().time().deltaTValue();
                                //checking if eddies convected outside of the box
                                ppr = (rotate1 & (pp_[i]-p0)) + p0;	//rotated to imaginary bounding box
                                if ( ppr[0] - sqrt(pow(sigma_[i][0],2)+pow(sigma_[i][1],2)) > endPosition[0] )
                                {
                        		Info << "generating new eddy"  << endl;
			                pp_[i] = ranGen_.position(startPosition,endPosition);	//generate point in fake coordinate
					//Info << "positions = " << pp_[i] << endl;
					pp_[i][0]=startPosition[0];
                                        sigma_[i] = ranGen_.position(sigmaMin_,sigmaMax_);
					intensity_[i] = fstu;
					if (pp_[i][1]<=delta_) ////
					{
                	                	intensity_[i] = vector(0,0,0) ;
					}

                 			pp_[i] = (rotate2 & (pp_[i]-p0)) + p0; //rotate it to the originat coordinate
                              		ranGen_.randomise(tempsign);
                                	tempsign = tempsign - tempvector;
                                	rndsign_[i][0]=sign(tempsign[0]);
                                	rndsign_[i][1]=sign(tempsign[1]);
                                	rndsign_[i][2]=sign(tempsign[2]);


		
                                }
			

                        }



                        pointsDict_.set("points",pp_);
                        pointsDict_.set("sigma",sigma_);
                        pointsDict_.set("intensity",intensity_);
                        pointsDict_.set("rndsign",rndsign_);

                        pointsDict_.set("time",this->db().time().timeOutputValue());

			//Info << "positions = " << pp_ << endl;
                        //Info << "sigma = " << sigma_ << endl;
                        //Info << "intensity = " << intensity_ << endl;
                        //Info << "rndsign = " << rndsign_ << endl;

                        //pointsDict_.Foam::regIOobject::write();
                        //Info << "points = " << pp_ << endl;
                }



                //Info << "out of else " << endl;
                //Info << "positions = " << pp << endl;
                //Info << "sigma = " << sigma << endl;
                //Info << "intensity = " << intensity << endl;
                //Info << "rndsign = " << rndsign << endl;


                scalar Vb((max(pp_ + sigma_)[0]-min(pp_ - sigma_)[0])*(endPosition[1]-bb.min()[1])*(bb.max()[2]-bb.min()[2]));
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

	        vectorField Uinlet = this->patch().Cf() ;

                forAll ( c,facei )
                {
                        patchField[facei] = pTraits<Type>::zero;//vector(0,0,0);
                        forAll(pp_,i)
                        {
                                dx = c[facei]-pp_[i];
                                dx[0]=dx[0]/sigma_[i][0];
                                dx[1]=dx[1]/sigma_[i][1];
                                dx[2]=dx[2]/sigma_[i][2];


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
                                        signedintensity[i][j] = rndsign_[i][j] * intensity_[i][j];
                                }
                                //patchField[facei] = unit;
				
				
				signedintensityType[i] = simplify(signedintensity[i], pTraits<Type>::one);
				patchField[facei] = patchField[facei] +
                                                                        (
                                                                        (signedintensityType[i])*mag(Uinf_)
                                                                        /sqrt(sigma_[i][0]*sigma_[i][1]*sigma_[i][2])
                                                                        * f
                                                                        * sqrt(Vb) / sqrt(n_)
                                                                        );

                                //scalar bb(0);
                                //bb = vector(sigma[i][0],0,0) & (vector(0,sigma[i][1],0) ^ vector(0,0,sigma[i][2]));//sigma[i] & vector(1,1,1);//a & signedintensity[i];
                                //Info << "result = " << sigma[i] << endl; //   /(sigma(1)*sigma(2)*sigma(3))* sqrt(Vb) *f((evalPosition - eddyPosition)./sigma);       

                        }





/*
valueExpression "d < 0.002 ? vector(1.785*d/0.002,5.635*d/0.002,0):vector(780.95*pow(d,3) - 431.07*pow(d,2) + 91.56*d + 1.6062,-175.24*pow(d,3) + 108.65*pow(d,2) - 32.449*d + 5.7001,0)";
        variables       "xp=0.0075005191;yp=0.0546987521;d=pow( pow(pos().x-xp,2)+ pow(pos().y-yp,2) ,0.5);";
*/
			xf = c[facei][0]; // face i, x location
			yf = c[facei][1]; // face i, y location

			d = sqrt( pow(xf-xp,2)+ pow(yf-yp,2) );		   		
			//interpolation on EXP profile
 			//Coeff[facei]=1;
			if ( d < 0.002 )
			{
				Uinlet[facei][0]=1.785*d/0.002;
				Uinlet[facei][1]=5.635*d/0.002;
				Uinlet[facei][2]=0;

			}
		 	else
			{
				Uinlet[facei][0]=780.95*pow(d,3) - 431.07*pow(d,2) + 91.56*d + 1.6062;
                                Uinlet[facei][1]=-175.24*pow(d,3) + 108.65*pow(d,2) - 32.449*d + 5.7001;
                                Uinlet[facei][2]=0;

			}	

			patchField[facei] = patchField[facei] + simplify(Uinlet[facei],pTraits<Type>::one);


                        //Info << "c[facei] =  " << c[facei] << endl;
                        //Info << "facei    =  " << facei << endl;
                        //Info << "result2  =  " << patchField[facei] << endl;
                }
		Info << "Done!" << endl;
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
void turbulentInflowAmesFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("n") << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaMin") << sigmaMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaMax") << sigmaMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("intensityMin") << intensityMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("intensityMax") << intensityMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("Uinf") << Uinf_ << token::END_STATEMENT << nl;
    os.writeKeyword("L") << L_ << token::END_STATEMENT << nl;
    os.writeKeyword("delta") << delta_ << token::END_STATEMENT << nl;
    os.writeKeyword("maxy") << maxy_ << token::END_STATEMENT << nl;
    os.writeKeyword("fstu") << fstu_ << token::END_STATEMENT << nl;
    referenceField_.writeEntry("referenceField", os);
    pointsDict_.Foam::regIOobject::write();
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

																																		// ************************************************************************* //
