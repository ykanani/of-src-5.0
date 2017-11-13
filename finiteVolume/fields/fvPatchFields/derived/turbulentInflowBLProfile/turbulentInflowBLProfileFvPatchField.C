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

#include "turbulentInflowBLProfileFvPatchField.H"
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
turbulentInflowBLProfileFvPatchField<Type>::turbulentInflowBLProfileFvPatchField
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
        intensityWall_(1,1,1),	   	//max intensity in all direction                
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
turbulentInflowBLProfileFvPatchField<Type>::turbulentInflowBLProfileFvPatchField
(
    const turbulentInflowBLProfileFvPatchField<Type>& ptf,
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
        intensityWall_(ptf.intensityWall_),
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
turbulentInflowBLProfileFvPatchField<Type>::turbulentInflowBLProfileFvPatchField
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
        intensityWall_(dict.lookup("intensityWall")),
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
turbulentInflowBLProfileFvPatchField<Type>::turbulentInflowBLProfileFvPatchField
(
    const turbulentInflowBLProfileFvPatchField<Type>& ptf
)
:
    	fixedValueFvPatchField<Type>(ptf),
	pointsDict_(ptf.pointsDict_),
	ranGen_(ptf.ranGen_),
        n_(ptf.n_),
        sigmaMin_(ptf.sigmaMin_),
        sigmaMax_(ptf.sigmaMax_),
        intensityWall_(ptf.intensityWall_),
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
turbulentInflowBLProfileFvPatchField<Type>::turbulentInflowBLProfileFvPatchField
(
    const turbulentInflowBLProfileFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
	pointsDict_(ptf.pointsDict_),
	ranGen_(ptf.ranGen_),
        n_(ptf.n_),
        sigmaMin_(ptf.sigmaMin_),
        sigmaMax_(ptf.sigmaMax_),
        intensityWall_(ptf.intensityWall_),
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
void turbulentInflowBLProfileFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    referenceField_.autoMap(m);
}


template<class Type>
void turbulentInflowBLProfileFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const turbulentInflowBLProfileFvPatchField<Type>& tiptf =
        refCast<const turbulentInflowBLProfileFvPatchField<Type> >(ptf);

    referenceField_.rmap(tiptf.referenceField_, addr);
}


template<class Type>
void turbulentInflowBLProfileFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

	const objectRegistry& obr = this->db();

        const IOdictionary& transportProperties= obr.lookupObject<IOdictionary>("transportProperties");
	dimensionedScalar nu(transportProperties.lookup("nu"));	
	Info << "nu =" << nu <<  endl;

	Field<Type>& patchField = *this;
	Info << "in update coef" << endl;
	//reading and updating eddy positions and calculating patch velocity if it needs an update
	if (curTimeIndex_ != this->db().time().timeIndex()      )
        {
		//experiments profile
		int n_inter;
		n_inter = 201;
		float ubla[201]={ 
					0	,
					0.01615	,
					0.0323	,
					0.04845	,
					0.064599674	,
					0.08074837	,
					0.09689511	,
					0.11303859	,
					0.129177182	,
					0.145308932	,
					0.161431564	,
					0.17754248	,
					0.193638762	,
					0.209717175	,
					0.22577417	,
					0.241805891	,
					0.257808173	,
					0.273776555	,
					0.289706281	,
					0.305592313	,
					0.321429333	,
					0.337211758	,
					0.352933746	,
					0.368589213	,
					0.384171841	,
					0.399675095	,
					0.415092235	,
					0.430416336	,
					0.445640305	,
					0.460756896	,
					0.475758733	,
					0.490638331	,
					0.505388116	,
					0.520000449	,
					0.53446765	,
					0.548782019	,
					0.562935868	,
					0.57692154	,
					0.59073144	,
					0.60435806	,
					0.617794004	,
					0.63103202	,
					0.644065023	,
					0.656886123	,
					0.669488652	,
					0.681866191	,
					0.694012592	,
					0.705922008	,
					0.717588911	,
					0.729008118	,
					0.740174812	,
					0.75108456	,
					0.761733333	,
					0.772117521	,
					0.782233949	,
					0.792079889	,
					0.80165307	,
					0.810951688	,
					0.819974411	,
					0.828720383	,
					0.837189227	,
					0.845381045	,
					0.853296412	,
					0.860936374	,
					0.868302441	,
					0.875396574	,
					0.882221179	,
					0.888779088	,
					0.895073547	,
					0.901108199	,
					0.906887062	,
					0.912414514	,
					0.917695266	,
					0.922734342	,
					0.927537055	,
					0.932108984	,
					0.936455944	,
					0.940583966	,
					0.944499267	,
					0.948208225	,
					0.951717354	,
					0.955033278	,
					0.958162704	,
					0.961112398	,
					0.963889161	,
					0.966499806	,
					0.968951135	,
					0.971249916	,
					0.973402866	,
					0.975416627	,
					0.977297755	,
					0.979052694	,
					0.980687772	,
					0.982209176	,
					0.983622948	,
					0.984934969	,
					0.986150951	,
					0.987276429	,
					0.988316754	,
					0.989277086	,
					0.99016239	,
					0.990977435	,
					0.991726789	,
					0.992414818	,
					0.993045689	,
					0.993623369	,
					0.994151625	,
					0.994634031	,
					0.995073966	,
					0.995474623	,
					0.995839011	,
					0.99616996	,
					0.996470126	,
					0.996741999	,
					0.996987908	,
					0.997210026	,
					0.997410377	,
					0.997590846	,
					0.997753179	,
					0.997898998	,
					0.9980298	,
					0.99814697	,
					0.99825178	,
					0.998345405	,
					0.998428921	,
					0.998503316	,
					0.998569493	,
					0.998628277	,
					0.99868042	,
					0.998726608	,
					0.998767463	,
					0.998803551	,
					0.998835381	,
					0.998863417	,
					0.998888076	,
					0.998909734	,
					0.99892873	,
					0.998945366	,
					0.998959915	,
					0.998972621	,
					0.998983701	,
					0.998993349	,
					0.99900174	,
					0.999009025	,
					0.999015342	,
					0.999020811	,
					0.999025539	,
					0.999029621	,
					0.99903314	,
					0.999036169	,
					0.999038773	,
					0.999041008	,
					0.999042924	,
					0.999044563	,
					0.999045964	,
					0.999047159	,
					0.999048178	,
					0.999049044	,
					0.999049781	,
					0.999050405	,
					0.999050934	,
					0.999051382	,
					0.99905176	,
					0.999052079	,
					0.999052347	,
					0.999052573	,
					0.999052762	,
					0.999052921	,
					0.999053054	,
					0.999053165	,
					0.999053258	,
					0.999053335	,
					0.9990534	,
					0.999053453	,
					0.999053497	,
					0.999053534	,
					0.999053565	,
					0.99905359	,
					0.99905361	,
					0.999053627	,
					0.999053641	,
					0.999053653	,
					0.999053662	,
					0.99905367	,
					0.999053676	,
					0.999053681	,
					0.999053686	,
					0.999053689	,
					0.999053692	,
					0.999053694	,
					0.999053696	,
					0.999053697	,
					0.999053698	,
					0.999053699	,
					0.9990537	,
					0.999053701	,
					0.999053701	,
					0.999053702	,
					0.999053702	,
					0.999053702	,
					0.999053702};
					
		float etaRef[201]={
					0	,
					0.05	,
					0.1	,
					0.15	,
					0.2	,
					0.25	,
					0.3	,
					0.35	,
					0.4	,
					0.45	,
					0.5	,
					0.55	,
					0.6	,
					0.65	,
					0.7	,
					0.75	,
					0.8	,
					0.85	,
					0.9	,
					0.95	,
					1	,
					1.05	,
					1.1	,
					1.15	,
					1.2	,
					1.25	,
					1.3	,
					1.35	,
					1.4	,
					1.45	,
					1.5	,
					1.55	,
					1.6	,
					1.65	,
					1.7	,
					1.75	,
					1.8	,
					1.85	,
					1.9	,
					1.95	,
					2	,
					2.05	,
					2.1	,
					2.15	,
					2.2	,
					2.25	,
					2.3	,
					2.35	,
					2.4	,
					2.45	,
					2.5	,
					2.55	,
					2.6	,
					2.65	,
					2.7	,
					2.75	,
					2.8	,
					2.85	,
					2.9	,
					2.95	,
					3	,
					3.05	,
					3.1	,
					3.15	,
					3.2	,
					3.25	,
					3.3	,
					3.35	,
					3.4	,
					3.45	,
					3.5	,
					3.55	,
					3.6	,
					3.65	,
					3.7	,
					3.75	,
					3.8	,
					3.85	,
					3.9	,
					3.95	,
					4	,
					4.05	,
					4.1	,
					4.15	,
					4.2	,
					4.25	,
					4.3	,
					4.35	,
					4.4	,
					4.45	,
					4.5	,
					4.55	,
					4.6	,
					4.65	,
					4.7	,
					4.75	,
					4.8	,
					4.85	,
					4.9	,
					4.95	,
					5	,
					5.05	,
					5.1	,
					5.15	,
					5.2	,
					5.25	,
					5.3	,
					5.35	,
					5.4	,
					5.45	,
					5.5	,
					5.55	,
					5.6	,
					5.65	,
					5.7	,
					5.75	,
					5.8	,
					5.85	,
					5.9	,
					5.95	,
					6	,
					6.05	,
					6.1	,
					6.15	,
					6.2	,
					6.25	,
					6.3	,
					6.35	,
					6.4	,
					6.45	,
					6.5	,
					6.55	,
					6.6	,
					6.65	,
					6.7	,
					6.75	,
					6.8	,
					6.85	,
					6.9	,
					6.95	,
					7	,
					7.05	,
					7.1	,
					7.15	,
					7.2	,
					7.25	,
					7.3	,
					7.35	,
					7.4	,
					7.45	,
					7.5	,
					7.55	,
					7.6	,
					7.65	,
					7.7	,
					7.75	,
					7.8	,
					7.85	,
					7.9	,
					7.95	,
					8	,
					8.05	,
					8.1	,
					8.15	,
					8.2	,
					8.25	,
					8.3	,
					8.35	,
					8.4	,
					8.45	,
					8.5	,
					8.55	,
					8.6	,
					8.65	,
					8.7	,
					8.75	,
					8.8	,
					8.85	,
					8.9	,
					8.95	,
					9	,
					9.05	,
					9.1	,
					9.15	,
					9.2	,
					9.25	,
					9.3	,
					9.35	,
					9.4	,
					9.45	,
					9.5	,
					9.55	,
					9.6	,
					9.65	,
					9.7	,
					9.75	,
					9.8	,
					9.85	,
					9.9	,
					9.95	,
					10};
								
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
                                	intensity_[i] = tempi*0;//(intensityWall_- tempi)*(1-pp_[i][1]/delta_)+tempi;
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
                                	intensity_[i] = (intensityWall_- tempi)*(1-pp_[i][1]/delta_)+tempi;
				}

                               	
					//intensity_[i] = (intensityWall_-tempi)*(1-pp_[i][1]/0.0043)+tempi;
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
	
		scalar x0(1/(nu.value()*25));
		scalar eta(0);
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
				eta = c[facei][1]*sqrt(mag(Uinf_)/(nu.value()*x0)); 
				if (eta<=etaRef[k])
				{
					Coeff[facei]=((ubla[k]-ubla[k-1])/(etaRef[k]-etaRef[k-1])*(eta-etaRef[k-1])+ubla[k-1])*mag(Uinf_);
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
void turbulentInflowBLProfileFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("n") << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaMin") << sigmaMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaMax") << sigmaMax_ << token::END_STATEMENT << nl;
    os.writeKeyword("intensityWall") << intensityWall_ << token::END_STATEMENT << nl;
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
