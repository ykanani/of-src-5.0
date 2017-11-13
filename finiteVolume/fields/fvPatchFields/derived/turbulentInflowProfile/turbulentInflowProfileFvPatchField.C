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
		//experiments profile
		int n_inter;
		n_inter = 36;
		float uexp[36]={ 0,
					0.5066411,
					0.52205956,
					0.53850603,
					0.5588069,
					0.57782364,
					0.59684145,
					0.61277515,
					0.6305086,
					0.6464436,
					0.6675202,
					0.6847422,
					0.7019657,
					0.7222743,
					0.7366708,
					0.7528674,
					0.7718927,
					0.7860334,
					0.7996612,
					0.8135454,
					0.8279448,
					0.8459438,
					0.85725814,
					0.8701171,
					0.8834907,
					0.89763576,
					0.9130669,
					0.92361236,
					0.93595916,
					0.9472784,
					0.95834213,
					0.9701783,
					0.98356396,
					0.9905184,
					0.99670696,
					1};

		float yexp[36]={		0,
						4.94747E-05,
						6.37923E-05,
						8.17234E-05,
						9.96079E-05,
						0.000128386,
						0.000175294,
						0.000211361,
						0.000265536,
						0.000323359,
						0.000421006,
						0.000515073,
						0.000634522,
						0.000779316,
						0.000891547,
						0.001025512,
						0.001195703,
						0.001326068,
						0.001474568,
						0.001612187,
						0.001771556,
						0.001967141,
						0.002101166,
						0.002282309,
						0.002474324,
						0.002677208,
						0.002901832,
						0.003068499,
						0.003278657,
						0.003492453,
						0.003731633,
						0.003996187,
						0.00438763,
						0.004659495,
						0.005014766,
						0.005355573};

/*	float yexp[36]={			0,
				2.47374E-05,
			3.18962E-05,
			4.08617E-05,
			4.9804E-05,
			0.000064193,
			0.000087647,
			0.000105681,
			0.000132768,
			0.00016168,
			0.000210503,
			0.000257537,
			0.000317261,
			0.000389658,
			0.000445774,
			0.000512756,
			0.000597852,
			0.000663034,
			0.000737284,
			0.000806094,
			0.000885778,
			0.000983571,
			0.001050583,
			0.001141155,
			0.001237162,
			0.001338604,
			0.001450916,
			0.00153425,
			0.001639329,
			0.001746227,
			0.001865817,
			0.001998094,
			0.002193815,
			0.002329748,
			0.002507383,
			0.002677787};
*/

/*	float yexp[36]={ 	  0,
			1.23687E-05,
			1.59481E-05,
			2.04309E-05,
			0.000024902,
			3.20965E-05,
			4.38235E-05,
			5.28405E-05,
			0.000066384,
			0.00008084,
			0.000105252,
			0.000128769,
			0.000158631,
			0.000194829,
			0.000222887,
			0.000256378,
			0.000298926,
			0.000331517,
			0.000368642,
			0.000403047,
			0.000442889,
			0.000491786,
			0.000525292,
			0.000570578,
			0.000618581,
			0.000669302,
			0.000725458,
			0.000767125,
			0.000819665,
			0.000873114,
			0.000932909,
			0.000999047,
			0.001096908,
			0.001164874,
			0.001253692,
			0.001338894};
*/			
	// defining eddy data file
	
	// defining eddy data file 
                // getting bound box for patch
                boundBox bb(this->patch().patch().localPoints(), true);
                vector startPosition(bb.min()[0]-L_,bb.min()[1],bb.min()[2]);
//                vector endPosition(bb.min()[0],bb.max()[1],bb.max()[2]);
                vector endPosition(bb.min()[0],maxy_,bb.max()[2]);

                Info << "start = " << startPosition << endl;
                Info << "end = " << endPosition << endl;
                // constructing requierd variables
                //vectorField pp(n_);
                //vectorField sigma(n_);
                //vectorField intensity(n_);
                //vectorField rndsign(n_);
                vectorField signedintensity(n_);
		Field<Type> signedintensityType(n_);
		vector  tempvector(0.5,0.5,0.5);
                vector  tempsign(1,1,1);
                vector  unit(1,1,1);
		vector  tempi(fstu_,fstu_,fstu_);
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
                        forAll ( pp_,i )
                        {
                                pp_[i] = ranGen_.position(startPosition,endPosition);
                                
				//rndsign[i]=sign(tempsign);
                                sigma_[i] = ranGen_.position(sigmaMin_,sigmaMax_);
				intensity_[i] = tempi;
				if (pp_[i][1]<=delta_) ////
				{
                                	intensity_[i] = (intensityMax_- tempi)*(1-pp_[i][1]/delta_)+tempi;
				}
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
                        //reading positions             
                        //pp              = vectorList(pointsDict_.lookup("points"));
                        //sigma           =vectorList(pointsDict_.lookup("sigma"));
                        //intensity       =vectorList(pointsDict_.lookup("intensity"));
                        //rndsign         =vectorList(pointsDict_.lookup("rndsign"));

                        //convecting eddies with mean velocity
                        vector a(this->db().time().timeOutputValue(),1,1);   //dummy
                        forAll ( pp_,i )
                        {
                                pp_[i]=pp_[i]+Uinf_*this->db().time().deltaTValue();
                                //checking if eddies convected outside of the box
                                if ( pp_[i][0] - sigma_[i][0] > bb.min()[0] )
                                {
                                        //generating new eddy
                                        //Info << "generating new eddy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1" << endl;
                                        pp_[i] = ranGen_.position(startPosition,endPosition);
					pp_[i][0]=startPosition[0];
                                        sigma_[i] = ranGen_.position(sigmaMin_,sigmaMax_);
					intensity_[i] = tempi;
				if (pp_[i][1]<=delta_) ////
				{
                                	intensity_[i] = (intensityMax_- tempi)*(1-pp_[i][1]/delta_)+tempi;
				}

                               	
					//intensity_[i] = (intensityMax_-tempi)*(1-pp_[i][1]/0.0043)+tempi;
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

	        scalarField Coeff = this->patch().Cf() & m;

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
                        //adding mean velocity
                        /*
			scalar yOverDelta = ( c[facei][1]  )/0.0043 ;

                                    if (yOverDelta > 1.0)
                                            {
                                               Coeff[facei]= 1;

                                            }
                                    else
                                            {

                                                Coeff[facei]= 0.757* pow((c[facei][1]/0.0007), (1.0/6.7));

                                                }
 			*/
 			//interpolation on EXP profile
 			Coeff[facei]=1;
 			for (int k=1; k<n_inter; k=k+1)
			{
				if (c[facei][1]<=yexp[k])
				{
					Coeff[facei]=(uexp[k]-uexp[k-1])/(yexp[k]-yexp[k-1])*(c[facei][1]-yexp[k-1])+uexp[k-1];
					break;
				}
			}	

			patchField[facei] = patchField[facei] + simplify(Coeff[facei]*Uinf_,pTraits<Type>::one);


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
