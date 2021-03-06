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

#include "sem2FvPatchField.H"
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
sem2FvPatchField<Type>::sem2FvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    	fixedValueFvPatchField<Type>(p, iF),
	mapperPtr_(NULL),
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
        statisticsDict_
                        (
                            IOobject
                            (
                                "statistics",
                                this->db().time().constant(),//timeName(this->db().time().timeOutputValue()-this->db().time().deltaTValue()),
                                this->db(),
                                IOobject::MUST_READ,
                                IOobject::AUTO_WRITE
                            )
                        ),


	ranGen_(1),//clock().getTime()),
        n_(100),                        		//number of eddies
        sigmaMin_(0.5,0.5,0.5),         		//min length scale in all direction
        sigmaMax_(0.5,0.5,0.5),     	    	//max length scale in all direction     
        Uinf_(1,0,0),   			//mean velocity
	cf_(1),					//skin friction coeff
        L_(1),          				//eddy box length in flow direction
	referenceField_(p.size()),			//reference filed
        curTimeIndex_(-1),       			// time index 
	pp_(100),
        sigma_(100),
        intensity_(100),
        rndsign_(100),
	delta_(0.001),
	fstu_(0),
	maxy_(0.01),
	procEddyFace(Pstream::nProcs()),
	eddyFaceAssigned(0)

 

{}


template<class Type>
sem2FvPatchField<Type>::sem2FvPatchField
(
    const sem2FvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
	fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
	mapperPtr_(NULL),
        pointsDict_(ptf.pointsDict_),
	statisticsDict_(ptf.pointsDict_),
	ranGen_(1),//clock().getTime()),
        n_(ptf.n_),
        sigmaMin_(ptf.sigmaMin_),
        sigmaMax_(ptf.sigmaMax_),
        Uinf_(ptf.Uinf_),
        cf_(ptf.cf_),
        L_(ptf.L_),
	referenceField_(ptf.referenceField_, mapper),
        curTimeIndex_(-1),
        pp_(ptf.pp_),
        sigma_(ptf.sigma_),
        intensity_(ptf.intensity_),
        rndsign_(ptf.rndsign_),
	delta_(ptf.delta_),
	fstu_(ptf.fstu_),
	maxy_(ptf.maxy_),
	procEddyFace(ptf.procEddyFace),
	eddyFaceAssigned(ptf.eddyFaceAssigned)
{}


template<class Type>
sem2FvPatchField<Type>::sem2FvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    	fixedValueFvPatchField<Type>(p, iF),
        mapperPtr_(NULL),
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
        statisticsDict_
                        (
                            IOobject
                            (
                                "statistics",
                                this->db().time().constant(),//timeName(this->db().time().timeOutputValue()-this->db().time().deltaTValue()),
                                this->db(),
                                IOobject::MUST_READ,
                                IOobject::AUTO_WRITE
                            )
                        ),


	ranGen_(1),//clock().getTime()),
        n_(readScalar(dict.lookup("n"))),
        sigmaMin_(dict.lookup("sigmaMin")),
        sigmaMax_(dict.lookup("sigmaMax")),
        Uinf_(dict.lookup("Uinf")),
        cf_(readScalar(dict.lookup("cf"))),
        L_(readScalar(dict.lookup("L"))),
	referenceField_("referenceField", dict, p.size()),
        curTimeIndex_(-1),
        pp_(vectorList(pointsDict_.lookup("points"))),
        sigma_(vectorList(pointsDict_.lookup("sigma"))),
        intensity_(vectorList(pointsDict_.lookup("intensity"))),
        rndsign_(vectorList(pointsDict_.lookup("rndsign"))),
	delta_(readScalar(dict.lookup("delta"))),
	fstu_(readScalar(dict.lookup("fstu"))),
 	maxy_(readScalar(dict.lookup("maxy"))),
	procEddyFace(Pstream::nProcs()),
	eddyFaceAssigned(0)

{
//	Info << "Clock.Time=" << clock().getTime() << endl;
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
sem2FvPatchField<Type>::sem2FvPatchField
(
    const sem2FvPatchField<Type>& ptf
)
:
    	fixedValueFvPatchField<Type>(ptf),
        mapperPtr_(NULL),
	pointsDict_(ptf.pointsDict_),
        statisticsDict_(ptf.pointsDict_),
	ranGen_(ptf.ranGen_),
        n_(ptf.n_),
        sigmaMin_(ptf.sigmaMin_),
        sigmaMax_(ptf.sigmaMax_),
        Uinf_(ptf.Uinf_),
        cf_(ptf.cf_),
        L_(ptf.L_),
	referenceField_(ptf.referenceField_),
        curTimeIndex_(-1),
        pp_(ptf.pp_),
        sigma_(ptf.sigma_),
        intensity_(ptf.intensity_),
        rndsign_(ptf.rndsign_),
	delta_(ptf.delta_),
	fstu_(ptf.fstu_),
 	maxy_(ptf.maxy_),
	procEddyFace(ptf.procEddyFace),
	eddyFaceAssigned(ptf.eddyFaceAssigned)


{}


template<class Type>
sem2FvPatchField<Type>::sem2FvPatchField
(
    const sem2FvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
        mapperPtr_(NULL),
	pointsDict_(ptf.pointsDict_),
        statisticsDict_(ptf.pointsDict_),
	ranGen_(ptf.ranGen_),
        n_(ptf.n_),
        sigmaMin_(ptf.sigmaMin_),
        sigmaMax_(ptf.sigmaMax_),
        Uinf_(ptf.Uinf_),
        cf_(ptf.cf_),
        L_(ptf.L_),
        referenceField_(ptf.referenceField_),
        curTimeIndex_(-1),
        pp_(ptf.pp_),
        sigma_(ptf.sigma_),
        intensity_(ptf.intensity_),
        rndsign_(ptf.rndsign_),
	delta_(ptf.delta_),
	fstu_(ptf.fstu_),
	maxy_(ptf.maxy_),
	procEddyFace(ptf.procEddyFace),
	eddyFaceAssigned(ptf.eddyFaceAssigned)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void sem2FvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    referenceField_.autoMap(m);
    // Clear interpolator
    mapperPtr_.clear();
}


template<class Type>
void sem2FvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const sem2FvPatchField<Type>& tiptf =
        refCast<const sem2FvPatchField<Type> >(ptf);

    referenceField_.rmap(tiptf.referenceField_, addr);
    // Clear interpolator
    mapperPtr_.clear();
}


template<class Type>
void sem2FvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
	
	Field<Type>& patchField = *this;
	
	//Sout << "p." << Pstream::myProcNo() <<  "VelocityField = " << patchField << endl;
	//const polyMesh& pMesh = this->patch().boundaryMesh().mesh();

	//Info << "in update coef" << endl;

	bool verbos(false);
	scalar utau(sqrt(cf_*0.5*Uinf_[0]*Uinf_[0]));
	if (curTimeIndex_ == -1 )
	{
                
                //Experiments profile
		Info << "Reading Input Statistics" << endl;  
		scalarField samplePoints(List<scalar>(statisticsDict_.lookup("points")));
		symmTensorField sampleRplus(List<symmTensor>(statisticsDict_.lookup("R")));
		scalarField sampleUplus(List<scalar>(statisticsDict_.lookup("U")));
		tensorField Lund(sampleRplus.size(), pTraits<tensor>::zero);
		
		scalarField sampleU(sampleUplus*utau);
		symmTensorField sampleR(sampleRplus*utau*utau);
                scalarField sampleY(samplePoints*delta_);
		//maxy_= 1 //maxy is already read from dictionary

		Info << "Uinf = " << Uinf_ << endl;
		Info << "delta = " << delta_ << endl;
		Info << "Friction velocity = " << utau << endl;
		Info << "maximum y =" << maxy_ << endl;
	        //Info << "U=" << sampleUplus << endl; 
		//Info << "R=" << sampleRplus << endl;
 	        //Info << "Points=" << samplePoints << endl; 

 
	}
       	//reading and updating eddy positions and calculating patch velocity if it needs an update
	if (curTimeIndex_ != this->db().time().timeIndex()      )
        {
			//initializing patchfield values to zero
		patchField = patchField * 0.0;

		const vectorField& c = this->patch().Cf();
		//Sout << "p." << Pstream::myProcNo() << "face size = " << c.size() << endl;
		//forAll ( c,facei )
		//{
		//	Sout << "face[" << facei << "]=" << c[facei] << endl;
		//}

		//experiments profile
		scalarField samplePoints(List<scalar>(statisticsDict_.lookup("points")));
		symmTensorField sampleRplus(List<symmTensor>(statisticsDict_.lookup("R")));
		scalarField sampleUplus(List<scalar>(statisticsDict_.lookup("U")));
		tensorField Lund(sampleRplus.size(), pTraits<tensor>::zero);
		
		scalarField sampleU(sampleUplus*utau);
		symmTensorField sampleR(sampleRplus*utau*utau);
                scalarField sampleY(samplePoints*delta_);

		Lund.replace(tensor::XX, sqrt(sampleR.component(symmTensor::XX)));
		Lund.replace(tensor::YX, sampleR.component(symmTensor::XY)/(Lund.component(tensor::XX)+ROOTVSMALL));
    		Lund.replace(tensor::ZX, sampleR.component(symmTensor::XZ)/(Lund.component(tensor::XX)+ROOTVSMALL));
		Lund.replace(tensor::YY, sqrt(sampleR.component(symmTensor::YY)-sqr(Lund.component(tensor::YX))));
    		Lund.replace(tensor::ZY, (sampleR.component(symmTensor::YZ) - Lund.component(tensor::YX)*Lund.component(tensor::ZX) )
					/(Lund.component(tensor::YY)+ROOTVSMALL));
		Lund.replace(tensor::ZZ, sqrt(sampleR.component(symmTensor::ZZ) - sqr(Lund.component(tensor::ZX))-sqr(Lund.component(tensor::ZY))));
		
		int n_inter;
		n_inter = sampleR.size();
                // getting bound box for patch
                boundBox bb(this->patch().patch().localPoints(), true);
                vector startPosition(bb.min()[0]-L_,bb.min()[1],bb.min()[2]);
                vector endPosition(bb.min()[0],min(maxy_,bb.max()[1]),bb.max()[2]);

                Info << "start = " << startPosition << endl;
                Info << "end = " << endPosition << endl;
                vectorField signedintensity(n_);
		
		vector  tempvector(0.5,0.5,0.5);
                vector  tempsign(1,1,1);
		vector  tempi(fstu_,fstu_,fstu_);
		//genrating random positions for the first time
                bool isFirst(readBool(pointsDict_.lookup("isFirst")));
		//List<List<List<scalar>>> procEddyFace(Pstream::nProcs());
			
		List<List<scalar>>& eddyFace= procEddyFace[Pstream::myProcNo()];
 		eddyFace.setSize(n_); 
                //Sout << "isFirst==" << isFirst << ",proc=" << Pstream::myProcNo() << endl;
		//Sout << "procEddyFace = " << procEddyFace << "proc=" << Pstream::myProcNo() << endl; 
		//Info << "isFirst is " << isFirst << endl;
		vector dx(0,0,0);  
                scalar f(0);
	        vector UConv(Uinf_);


                if (isFirst )
                {

                        Info << "Generating Eddies and adding faces to eddies..." << endl;
                        //generating eddy positions
                        vector a(this->db().time().timeOutputValue(),1,1);      //dummy
                        forAll ( pp_,i )
                        {
                                pp_[i] = ranGen_.position(startPosition,endPosition);
                                //Sout << "pp[" << i << "]=" << pp_[i] << ", proc=" << Pstream::myProcNo() << endl;
				//rndsign[i]=sign(tempsign);
                                sigma_[i] = ranGen_.position(sigmaMin_,sigmaMax_);
				intensity_[i] = tempi;
                                ranGen_.randomise(tempsign);
                                tempsign = tempsign - tempvector;
                                rndsign_[i][0]=sign(tempsign[0]);
                                rndsign_[i][1]=sign(tempsign[1]);
                                rndsign_[i][2]=sign(tempsign[2]);
				//assigning face indices to each eddy
				
				List<scalar>& eddyfacei = eddyFace[i];
				eddyFaceAssigned=true;

				forAll ( c,facei )
				{
					dx = c[facei]-pp_[i];

                                	dx[0]=dx[0]/sigma_[i][0];
                                	dx[1]=dx[1]/sigma_[i][1];
                                	dx[2]=dx[2]/sigma_[i][2];
					//Info << "p=" << i << ",face=" << facei << ",dx=" << dx << endl;
                                        if (mag(dx[1]) <= 1 &&  mag(dx[2]) <= 1 )
                                        {
						eddyfacei.append(facei);
                                        }
				}
		        }
		        //Gathering information from all processors
			//Pstream::gatherList(procEddyFace);
			//if (Pstream::master())
			//{
				//Info << "procEddyFace =" << procEddyFace ;
				//pointsDict_.set("eddyFace",procEddyFace);

			//}
			pointsDict_.set("points",pp_);
			pointsDict_.set("sigma",sigma_);
			pointsDict_.set("intensity",intensity_);
			pointsDict_.set("rndsign",rndsign_);
			pointsDict_.set("time",this->db().time().timeOutputValue());
			pointsDict_.set("isFirst",false);
			pointsDict_.Foam::regIOobject::write();
			
			Info << "First Points = " << pp_[0] << endl;
                }
		else
                {
			if (!eddyFaceAssigned)
			{
				Info << "Assigning eddies to faces..." << endl;
				forAll (pp_,i)	
				{
					List<scalar>& eddyfacei = eddyFace[i];
					forAll ( c,facei )
					{
						dx = c[facei]-pp_[i];

						dx[0]=dx[0]/sigma_[i][0];
						dx[1]=dx[1]/sigma_[i][1];
						dx[2]=dx[2]/sigma_[i][2];
						//Info << "p=" << i << ",face=" << facei << ",dx=" << dx << endl;
						if (mag(dx[1]) <= 1 &&  mag(dx[2]) <= 1 )
						{
							eddyfacei.append(facei);
						}
					}
				}
				eddyFaceAssigned=true;

			}
                        Info << "Conveting Eddies..." << endl;
                        //convecting eddies with mean velocity
                        vector a(this->db().time().timeOutputValue(),1,1);   //dummy
                        int counter(0);
			int counterConvect(0);
			bool eddyFaceneedUpdate(0);
			forAll ( pp_,i )
                        {
			/*	//interpolation on EXP profile to calc convection velocity
				UConv=Uinf_;
				for (int k=1; k<n_inter; k=k+1)
				{
					if (pp_[i][1]<=sampleY[k])
					{
						UConv[0]=(sampleU[k]-sampleU[k-1])/(sampleY[k]-sampleY[k-1])*(pp_[i][1]-sampleY[k-1])+sampleU[k-1];
						break;
					}
				}
			*/
				pp_[i]=pp_[i]+Uinf_*this->db().time().deltaTValue();
                                 //checking if eddies convected outside of the box
	                        if ( pp_[i][0] - sigma_[i][0] > bb.min()[0] )
                                {
                                	counterConvect++;
				        //generating new eddy
                                        //Info << "generating new eddy!!!!!!!!!!!!!!!!!1" << endl;
                                        pp_[i] = ranGen_.position(startPosition,endPosition);
					pp_[i][0]=startPosition[0];
                                        sigma_[i] = ranGen_.position(sigmaMin_,sigmaMax_);
					intensity_[i] = tempi;
                               	
					//intensity_[i] = (intensityMax_-tempi)*(1-pp_[i][1]/0.0043)+tempi;
                              		ranGen_.randomise(tempsign);
                                	tempsign = tempsign - tempvector;
                                	rndsign_[i][0]=sign(tempsign[0]);
                                	rndsign_[i][1]=sign(tempsign[1]);
                                	rndsign_[i][2]=sign(tempsign[2]);
					//assigning face indices to each eddy
					List<scalar>& eddyfacei = eddyFace[i];

					eddyfacei.clear();
					eddyFaceneedUpdate=true;
					forAll ( c,facei )
					{
						dx = c[facei]-pp_[i];

						dx[0]=dx[0]/sigma_[i][0];
						dx[1]=dx[1]/sigma_[i][1];
						dx[2]=dx[2]/sigma_[i][2];
						//Info << "p=" << i << ",face=" << facei << ",dx=" << dx << endl;

						if (mag(dx[1]) <= 1 &&  mag(dx[2]) <= 1 )
						{
							eddyfacei.append(facei);
						}
					}
					//Sout << "ppNo=" << i << "eddfacei.size=" << eddyfacei.size() << ",proc=" << Pstream::myProcNo() << endl;
						
                                }

		          }
			//if (eddyFaceneedUpdate)
			//{
			//	Pstream::gatherList(procEddyFace);
			//	//Pstream::scatterList(procEddyFace);
			//	if (Pstream::master())
			//	{
			//	//Info << "procEddyFace =" << procEddyFace ;
			//	//pointsDict_.set("eddyFace",procEddyFace);
			//	}
			//}

			pointsDict_.set("points",pp_);
			pointsDict_.set("sigma",sigma_);
			pointsDict_.set("intensity",intensity_);
			pointsDict_.set("rndsign",rndsign_);
			pointsDict_.set("time",this->db().time().timeOutputValue());
			//pointsDict_.Foam::regIOobject::write();
		
                        Info << "First Points = " << pp_[0] << endl;
                } // end of else


                        //if(verbos)
                        //{

			//	Info << "out of else " << endl;
			//	Info << "positions = " << pp << endl;
			//	Info << "sigma = " << sigma << endl;
			//	Info << "intensity = " << intensity << endl;
			//	Info << "rndsign = " << rndsign << endl;
			//}

                scalar Vb((max(pp_ + sigma_)[0]-min(pp_ - sigma_)[0])*(endPosition[1]-bb.min()[1])*(bb.max()[2]-bb.min()[2]));
                //Info << "maximum x = " << max(pp + sigma)[0] << endl;         
                //Info << "minimum x = " << min(pp - sigma)[0] << endl; 
                //Info << "Vb      = " << Vb << endl;
                tensor a(1,0,0,0,1,0,0,0,1);

                //calculating patch velocity based on eddies current poistion
                //vectorField& patchField = *this;
                Type test;
                vector Lcy(0,bb.max()[1]-bb.min()[1],0); //cyclic length in y
                vector Lcz(0,0,bb.max()[2]-bb.min()[2]); //cyclic length in z

		vector unit(1,1,1);
		vector unitx(1,0,0);
		Type unitType(pTraits<Type>::one);
	        vector m(0,0,0);
		vector ppPlusSigma(0,0,0);
                vector ppMinusSigma(0,0,0);
                vector Dbb(0,0,0); // Distance from bounding box
	        scalar UCoeff(1);
		tensor LundCoeff(pTraits<tensor>::zero);
		//Info << "starting eddy loop " << endl;
                forAll ( pp_,i )
                {
                        //patchField[facei] = pTraits<Type>::zero;//vector(0,0,0);
			//Info << "Eddy No. " << i << endl;
			
                        for(int jj=0; jj<eddyFace[i].size();jj=jj+1)
                        {
				int facei = eddyFace[i][jj];
                                dx = c[facei]-pp_[i];
                                dx[0]=dx[0]/sigma_[i][0];
                                dx[1]=dx[1]/sigma_[i][1];
                                dx[2]=dx[2]/sigma_[i][2];
				//Sout <<  "proc," << Pstream::myProcNo() << ",ppNo,"<< i << ",pp," << pp_[i] << ",facei," << facei << ",c," << c[facei] << ",dx," << dx << endl;
                                f=1;
                                for(int j=0; j<3; j=j+1)
                                {
                                        if (mag(dx[j]) >= 1 ) // no need for this, but just in case!
                                        {
                                                f = 0;
                                        }
                                        else
                                        {
                                                f = f * sqrt(scalar(1.5))*(1-mag(dx[j]));      //T21
                                        }
                                }
				
				patchField[facei] = patchField[facei] +
                                                                        (
                                                                        simplify(rndsign_[i], pTraits<Type>::one)
									/sqrt(sigma_[i][0]*sigma_[i][1]*sigma_[i][2])
                                                                        * f
                                                                        * sqrt(Vb) / sqrt(n_)
                                                                        );


 				//Ensuring cyclic behaviour
                                ppPlusSigma= pp_[i]+sigma_[i];
                                ppMinusSigma= pp_[i]-sigma_[i];
                                //check in z direction
                                if (ppPlusSigma[2] > bb.max()[2]) //outside bb from right
                                {
                                        dx = c[facei]-(pp_[i]-Lcz);
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
                                                        f = f * sqrt(scalar(1.5))*(1-mag(dx[j]));      //T21
                                                }
                                        }
					patchField[facei] = patchField[facei] +
                                                                        (
                                                                        simplify(rndsign_[i], pTraits<Type>::one)
                                                                        /sqrt(sigma_[i][0]*sigma_[i][1]*sigma_[i][2])
                                                                        * f
                                                                        * sqrt(Vb) / sqrt(n_)
                                                                        );


                                } //end if outside bb from right
                                else if (ppMinusSigma[2] < bb.min()[2]) //outside bb from left
                                {
                                        dx = c[facei]-(pp_[i]+Lcz);
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
                                                        f = f * sqrt(scalar(1.5))*(1-mag(dx[j]));      //T21
                                                }
                                        }
					patchField[facei] = patchField[facei] +
                                                                        (
                                                                        simplify(rndsign_[i], pTraits<Type>::one)
                                                                        /sqrt(sigma_[i][0]*sigma_[i][1]*sigma_[i][2])
                                                                        * f
                                                                        * sqrt(Vb) / sqrt(n_)
                                                                        );

                                }

                                 
                        }
		}
		Info << " Constructing mean field " << endl;	
		forAll(c,facei)
		{
                        //interpolation on EXP profile
 			UCoeff=Uinf_[0];
			LundCoeff=pTraits<tensor>::zero;
 			for (int k=1; k<n_inter; k=k+1)
			{
				if (c[facei][1]<=sampleY[k])
				{
					UCoeff=(sampleU[k]-sampleU[k-1])/(sampleY[k]-sampleY[k-1])*(c[facei][1]-sampleY[k-1])+sampleU[k-1];
					LundCoeff=(Lund[k]-Lund[k-1])/(sampleY[k]-sampleY[k-1])*(c[facei][1]-sampleY[k-1])+Lund[k-1];
					break;
				}
			}	

			patchField[facei] = dotp(LundCoeff,patchField[facei]) + simplify(UCoeff*unitx,pTraits<Type>::one);


                        //Info << "y/delta =  " << c[facei][1]/delta_ << endl;
                        //Info << "facei    =  " << facei << endl;
                        //Info << "UCoef  =  " << UCoeff << endl;
                        //Info << "LundCoef  =  " << LundCoeff << endl;
			
                }
		//Sout << "p=" << Pstream::myProcNo() << ", velocity = " << patchField << endl;
                curTimeIndex_ = this->db().time().timeIndex();
		Info << "Done!" << endl;

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
void sem2FvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("n") << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaMin") << sigmaMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaMax") << sigmaMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("Uinf") << Uinf_ << token::END_STATEMENT << nl;
    os.writeKeyword("cf") << cf_ << token::END_STATEMENT << nl;
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
