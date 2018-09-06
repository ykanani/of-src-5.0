/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::porosityModels::tubeBank::apply
(
    scalarField& Udiag,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const scalar sl = sl_;
    const scalar lp = lp_;
    const scalar st = st_;
    const scalar Dh = Dh_;
    const scalar d = d_;
	const scalar c0((1/sl)*0.5*0.938*pow(d/lp,-0.283)*pow((sl-lp)/Dh,0.171)*pow((st-d)/Dh,-1.19)*pow(Dh/1.5e-5,-0.274)*pow(st/(st-d),1.726));

    Info << "C0 is :" << c0 << endl;
    //const scalar Dh=pow(3.14*d_,2)+(4*(Lp_-d_)*d_)/(3.14*d_)+2*(Lp_-d_);
    //const scalar a=pow((sl_-Lp_)/Dh,0.134)*pow((st_-D_)/Dh,-1.19);
    //const scalar b=pow((st_-D_)/Dh,-1.19);
    //const scalar c=0.983*pow((D_/Lp_),-0.283);
    //const scalar e=pow((0.0254*st_)/((0.0254*st_)-D_),-0.274);
    //const scalar C0=C0_;

    forAll(cellZoneIDs_, zoneI)
    {
        
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];
             
        forAll(cells, i)
        {
            const label celli = cells[i];
            
            Udiag[celli] +=
                V[celli]*rho[celli]*c0*pow(magSqr(U[celli]), 0.363);
             
        }
      
	}
	
}


template<class RhoFieldType>
void Foam::porosityModels::tubeBank::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const scalar sl = sl_;
    const scalar st = st_;
    const scalar lp = lp_;
    const scalar Dh = Dh_;
    const scalar d = d_;
    const scalar c0((0.5*0.938*pow(d/lp,-0.283)*pow((sl-lp)/Dh,0.171)*pow((st-d)/Dh,-1.19)*pow(Dh/1.5e-5,-0.274)*pow(st/(st-d),1.726))/Dh);
    Info << "I AM BEING USED :" << c0 << endl;
    /*const scalar Dh=pow(3.14*d_,2)+(4*(Lp_-d_)*d_)/(3.14*d_)+2*(Lp_-d_);
    const scalar a=pow((sl_-Lp_)/Dh,0.134);
    const scalar b=pow((st_-D_)/Dh,-1.19);
    const scalar c=0.983*pow((D_/Lp_),-0.283);
    const scalar e=pow((0.0254*st_)/((0.0254*st_)-D_),-0.274);
    */

    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];
        
        forAll(cells, i)
        {
            const label celli = cells[i];

            AU[celli] =
                AU[celli] + I*(rho[celli]*c0*pow(magSqr(U[celli]), 0.363));
        }
    }
}


// ************************************************************************* //
