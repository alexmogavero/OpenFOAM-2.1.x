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

#include "exactRiemannSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//Equation 4.6 and 4.7
double exactRiemannSolver::f
(
	const double& p,
	const double& pLR,
	const double& rhoLR
) const
{
	double A(2/((gamma_+1)*rhoLR)); 		//eq. 4.8
	double B(((gamma_-1)/(gamma_+1))*pLR);  //eq. 4.8
	double aLR(::sqrt(gamma_*pLR/rhoLR));
	if(p>pLR) //Shock
	{
		return (p - pLR)*::sqrt(A/(p + B));
	}else  //Expansion
	{
		return (2*aLR/(gamma_-1))*(::pow((p/pLR),((gamma_-1)/(2*gamma_))) - 1);
	}
}

double exactRiemannSolver::fDer
(
	const double& p,
	const double& pLR,
	const double& rhoLR
) const
{
	double A(2/((gamma_+1)*rhoLR)); 		//eq. 4.8
	double B(((gamma_-1)/(gamma_+1))*pLR);  //eq. 4.8
	double aLR(::sqrt(gamma_*pLR/rhoLR));
	if(p>pLR) //Shock
	{
		return ::sqrt(A/(p + B))*(1 - (p - pLR)/(2*(B+p)));   //eq.4.37
	}else  //Expansion
	{
		return (::pow((p/pLR),(-gamma_-1)/(2*gamma_)) - 1)/(rhoLR*aLR);   //eq.4.37
	}
}

double
exactRiemannSolver::pStar
(
	double p
) const
{
	double val = 1000;
	double val00;
	double valDer;
	double p0 = 0.0;
	double p00;
	const double toll=1e-6;
	int n=0;
	const int nMax=50;
	double rP=1;
	double rP0=1;
	const double omega = 0.5;
	double rVal=1;

	if(exactRiemannSolver::debug)
	{
		Info << "exactRiemannSolver::pStar(double p)" << endl;
		Info << "p\t\trP\t\trVal" << endl;
		Info << p << " " << rP << " " << rVal << endl;
	}
	do
	{
		n++;
		val00 = val;
		val = f(p,pL_,rhoL_) + f(p,pR_,rhoR_) + uR_ - uL_; //eq. 4.5
		valDer = fDer(p,pL_,rhoL_) + fDer(p,pR_,rhoR_);
		p00 = p0;
		p0 = p;
		rP0 = rP;
		if(n<10)
		{
			p -= val/valDer;  //Newton Raphson
		}else
		{
			p = (val00*p0 - val*p00)/(val00 - val);  //Secant method
		}
		if(p<0)
		{
			p = toll;
		}
		rP = ::fabs(p-p0)/(0.5*(p+p0));
		rVal = 2*val/(uR_ + uL_ + toll);
		if(rP/rP0>0.1)
		{
			p = p0 + omega*(p - p0);  //underRelax
		}

		if(exactRiemannSolver::debug)
		{
			Info << p << " " << rP << " " << rVal << endl;
		}
	}while(rP>toll && n<nMax);

	if(n==nMax)
	{
		Info << "Warning: maximum number of iterations exceeded" << endl;
	}

    return p;
}

double
exactRiemannSolver::pStarTwoRar
(
) const
{
	double aL(::sqrt(gamma_*pL_/rhoL_));
	double aR(::sqrt(gamma_*pR_/rhoR_));
	double G((gamma_-1)/(2*gamma_));

	if(2*aL/(gamma_-1) + 2*aR/(gamma_-1)<=uR_ - uL_)   //eq.4.82
	{
		Info << "Warning: Vacuum creation" << endl;
		return 0.0;
	}else
	{
		return ::pow((aL + aR - 0.5*(gamma_-1)*(uR_ - uL_))/
				(aL/::pow(pL_,G) + aR/::pow(pR_,G)),1/G);		//Eq.4.46
	}
}

double
exactRiemannSolver::uStar
(
	double& p
) const
{
	return 0.5*(uL_ + uR_ + f(p,pR_,rhoR_) - f(p,pL_,rhoL_));  //eq. 4.9
}

void
exactRiemannSolver::uStar
(
)
{
	if(waveC_==VACUUM)
	{
		double aL(::sqrt(gamma_*pL_/rhoL_));
		double aR(::sqrt(gamma_*pR_/rhoR_));
		uStarL_ = uL_ + 2*aL/(gamma_-1);
		uStarR_ = uR_ - 2*aR/(gamma_-1);
	}else
	{
		uStarL_ = uStar(pStar_);
		uStarR_ = uStarL_;
	}
}

void exactRiemannSolver::setWave
(
	double& pStar
)
{
	if(pStar==0.0)
	{
		waveC_ = VACUUM;
	}else
	{
		waveC_ = CONTACT;
	}
	if(pStar>pL_)
	{
		waveL_ = SHOCK;
	}else
	{
		waveL_ = RAREFACTION;
	}

	if(pStar>pR_)
	{
		waveR_ = SHOCK;
	}else
	{
		waveR_ = RAREFACTION;
	}
}

void exactRiemannSolver::solve
(
)
{
	pStar_ = pStarTwoRar();
	setWave(pStar_);
	if(waveL_==SHOCK || waveR_==SHOCK)
	{
		pStar_ = pStar(pStar_);
		setWave(pStar_);
	}
	uStar();
	rhoStar();
	speed();

	if(exactRiemannSolver::debug==1)
	{
		Info << "exactRiemannSolver::solve()" << endl;
		Info << "waves: " << waveToString(waveL_) << " "
			 << waveToString(waveC_) << " " << waveToString(waveR_) << endl;
	}
}

void exactRiemannSolver::rhoStar
(
)
{
	double G = (gamma_ - 1)/(gamma_ + 1);

	if(waveC_==VACUUM)
	{
		rhoStarL_ = 0.0;
		rhoStarR_ = 0.0;
	}else
	{
		double pRatio(pStar_/pL_);
		if(waveL_==SHOCK)
		{
			rhoStarL_ = rhoL_*(pRatio + G)/(G*pRatio + 1); //eq.4.50
		}else if(waveL_==RAREFACTION)
		{
			rhoStarL_ = rhoL_*::pow(pRatio,1/gamma_);	//Eq.4.53
		}

		pRatio = pStar_/pR_;
		if(waveR_==SHOCK)
		{
			rhoStarR_ = rhoR_*(pRatio + G)/(G*pRatio + 1); //eq.4.57
		}else if(waveR_==RAREFACTION)
		{
			rhoStarR_ = rhoR_*::pow(pRatio,1/gamma_);	//eq.4.60
		}
	}
}

void exactRiemannSolver::speed
(
)
{
	double pRatio(pStar_/pL_);
	double aL(::sqrt(gamma_*pL_/rhoL_));
	if(waveL_==SHOCK)
	{
		speed_[0] = uL_ - aL*::sqrt(pRatio*(gamma_+1)/(2*gamma_) +
				(gamma_ - 1)/(2*gamma_));  //eq.4.52
		speed_[1] = speed_[0];
	}else if(waveL_==RAREFACTION)
	{
		double aStarL;
		if (waveC_==VACUUM)
		{
			aStarL = 0.0;
		}else
		{
			aStarL = ::sqrt(gamma_*pStar_/rhoStarL_);
		}
		speed_[0] = uL_ - aL;  //eq.4.55
		speed_[1] = uStarL_ - aStarL;   //eq.4.55
	}

	speed_[2] = 0.5*(uStarL_ + uStarR_);

	pRatio = pStar_/pR_;
	double aR(::sqrt(gamma_*pR_/rhoR_));
	if(waveR_==SHOCK)
	{
		speed_[3] = uR_ + aR*::sqrt(pRatio*(gamma_+1)/(2*gamma_) +
				(gamma_ - 1)/(2*gamma_));  //eq.4.59
		speed_[4] = speed_[3];
	}else if(waveR_==RAREFACTION)
	{
		double aStarR;
		if (waveC_==VACUUM)
		{
			aStarR = 0.0;
		}else
		{
			aStarR = ::sqrt(gamma_*pStar_/rhoStarR_);
		}
		speed_[3] = uStarR_ + aStarR;	//eq.4.62
		speed_[4] = uR_ + aR;			//eq.4.62
	}

	//TODO to decide if it is worth to keep this check, since it may be resources consuming
	double toll = 1e-6;
	for(int i=0; i<4; i++)
	{
		if(speed_[i+1]-speed_[i] < -toll*::fabs(speed_[i]))
		{
			Info << "Warning: Characteristic speeds not sorted." << endl;
			Info << speed_[0] << " " << speed_[1] << " " << speed_[2] << " " << speed_[3] << " " << speed_[4] << endl;
			break;
		}
	}
}

double exactRiemannSolver::p
(
	double xt
)const
{
	double p(0);

	if(xt<=speed_[0])
	{
		p = pL_;
	}else if(xt>speed_[0] && xt<speed_[1])
	{
		switch(waveL_)
		{
		case SHOCK:
			Info << "Warning: wave tagged as SHOCK has 2 different propagation speeds" << endl;
			p = pL_;
			break;
		case RAREFACTION:
		{
			double aL(::sqrt(gamma_*pL_/rhoL_));
			p = pL_*::pow((2/(gamma_+1) + (uL_-xt)*(gamma_-1)/((gamma_+1)*aL)),2*gamma_/(gamma_-1));   //eq.4.56
			break;
		}
		case NONE:
			FatalErrorIn
			(
				"exactRiemannSolver::p="
				"(double& xt)const"
			)   << "waveL must be assigned!"
				<< exit(FatalError);
				break;
		default:
			FatalErrorIn
				(
					"exactRiemannSolver::u="
					"(double& xt)const"
				)   << "waveL uncorrect type!"
					<< exit(FatalError);
			break;
		}
	}else if(xt>=speed_[1] && xt<=speed_[3])
	{
		p = pStar_;
	}else if(xt>speed_[3] && xt<speed_[4])
	{
		switch(waveR_)
		{
		case SHOCK:
			Info << "Warning: wave tagged as SHOCK has 2 different propagation speeds" << endl;
			p = pR_;
		break;
		case RAREFACTION:
		{
			double aR(::sqrt(gamma_*pR_/rhoR_));
			p = pR_*::pow((2/(gamma_+1) - (uR_-xt)*(gamma_-1)/((gamma_+1)*aR)),2*gamma_/(gamma_-1));   //eq.4.63
			break;
		}
		case NONE:
			FatalErrorIn
			(
				"exactRiemannSolver::p="
				"(double& xt)const"
			)   << "waveL must be assigned!"
				<< exit(FatalError);
				break;
		default:
			FatalErrorIn
				(
					"exactRiemannSolver::u="
					"(double& xt)const"
				)   << "waveL uncorrect type!"
					<< exit(FatalError);
			break;
		}
	}else if(xt>=speed_[4])
	{
		p = pR_;
	}

	return p;
}

double exactRiemannSolver::u
(
	double xt
)const
{
	double u(0);

	if(xt<=speed_[0])
	{
		u = uL_;
	}else if(xt>speed_[0] && xt<speed_[1])
	{
		switch(waveL_)
		{
		case SHOCK:
			Info << "Warning: wave tagged as SHOCK has 2 different propagation speeds" << endl;
			u = uL_;
			break;
		case RAREFACTION:
		{
			double aL(::sqrt(gamma_*pL_/rhoL_));
			u = (2/(gamma_+1))*(aL + 0.5*(gamma_-1)*uL_ + xt);  //eq.4.56
			break;
		}
		case NONE:
			FatalErrorIn
			(
				"exactRiemannSolver::u="
				"(double& xt)const"
			)   << "waveL must be assigned!"
				<< exit(FatalError);
				break;
		default:
			FatalErrorIn
				(
					"exactRiemannSolver::u="
					"(double& xt)const"
				)   << "waveL uncorrect type!"
					<< exit(FatalError);
			break;
		}
	}else if(xt>=speed_[1] && xt<speed_[2])
	{
		u = uStarL_;
	}else if(xt>=speed_[2] && xt<=speed_[3])
	{
		u = uStarR_;
	}else if(xt>speed_[3] && xt<speed_[4])
	{
		switch(waveR_)
		{
		case SHOCK:
			Info << "Warning: wave tagged as SHOCK has 2 different propagation speeds" << endl;
			u = uR_;
		break;
		case RAREFACTION:
		{
			double aR(::sqrt(gamma_*pR_/rhoR_));
			u = (2/(gamma_+1))*(-aR + 0.5*(gamma_-1)*uR_ + xt);  //eq.4.63
			break;
		}
		case NONE:
			FatalErrorIn
			(
				"exactRiemannSolver::u="
				"(double& xt)const"
			)   << "waveL must be assigned!"
				<< exit(FatalError);
				break;
		default:
			FatalErrorIn
				(
					"exactRiemannSolver::u="
					"(double& xt)const"
				)   << "waveL uncorrect type!"
					<< exit(FatalError);
			break;
		}
	}else if(xt>=speed_[4])
	{
		u = uR_;
	}

	return u;
}

double exactRiemannSolver::rho
(
	double xt
)const
{
	double rho(0);

	if(xt<=speed_[0])
	{
		rho = rhoL_;
	}else if(xt>speed_[0] && xt<speed_[1])
	{
		switch(waveL_)
		{
		case SHOCK:
			Info << "Warning: wave tagged as SHOCK has 2 different propagation speeds" << endl;
			rho = rhoL_;
			break;
		case RAREFACTION:
		{
			double aL(::sqrt(gamma_*pL_/rhoL_));
			rho = rhoL_*::pow((2/(gamma_+1) + (uL_-xt)*(gamma_-1)/((gamma_+1)*aL)),2/(gamma_-1));  //eq.4.56
			break;
		}
		case NONE:
			FatalErrorIn
			(
				"exactRiemannSolver::u="
				"(double& xt)const"
			)   << "waveL must be assigned!"
				<< exit(FatalError);
				break;
		default:
			FatalErrorIn
			(
				"exactRiemannSolver::u="
				"(double& xt)const"
			)   << "waveL uncorrect type!"
				<< exit(FatalError);
				break;
		}
	}else if(xt>=speed_[1] && xt<speed_[2])
	{
		rho = rhoStarL_;
	}else if(xt>=speed_[2] && xt<=speed_[3])
	{
		rho = rhoStarR_;
	}else if(xt>speed_[3] && xt<speed_[4])
	{
		switch(waveR_)
		{
		case SHOCK:
			Info << "Warning: wave tagged as SHOCK has 2 different propagation speeds" << endl;
			rho = rhoR_;
		break;
		case RAREFACTION:
		{
			double aR(::sqrt(gamma_*pR_/rhoR_));
			rho = rhoR_*::pow((2/(gamma_+1) - (uR_-xt)*(gamma_-1)/((gamma_+1)*aR)),2/(gamma_-1));  //eq.4.63
			break;
		}
		case NONE:
			FatalErrorIn
			(
				"exactRiemannSolver::u="
				"(double& xt)const"
			)   << "waveL must be assigned!"
				<< exit(FatalError);
		break;
		default:
			FatalErrorIn
			(
				"exactRiemannSolver::u="
				"(double& xt)const"
			)   << "waveL uncorrect type!"
				<< exit(FatalError);
		break;
		}
	}else if(xt>=speed_[4])
	{
		rho = rhoR_;
	}

	return rho;
}

double exactRiemannSolver::speed
(
	int i
)const
{
	return speed_[i];
}

void exactRiemannSolver::test
(
)
{
	const double rhoL[5] = {1.0,1.0,1.0,1.0,5.99924};
	const double uL[5] = {0.0,-2.0,0.0,0.0,19.5975};
	const double pL[5] = {1.0,0.4,1000.0,0.01,460.894};
	const double rhoR[5] = {0.125,1.0,1.0,1.0,5.99242};
	const double uR[5] = {0.0,2.0,0.0,0.0,-6.19633};
	const double pR[5] = {0.1,0.4,0.01,100.0,46.0950};
	for(int i=0;i<5;i++)
	{
		exactRiemannSolver reinmann(pL[i],pR[i],rhoL[i],rhoR[i],uL[i],uR[i],1.4);
		reinmann.solve();

		Info << i+1 << "--"
				<< " " << reinmann.p(0.5*(reinmann.speed(1)+reinmann.speed(2)))
				<< " " << reinmann.u(0.5*(reinmann.speed(1)+reinmann.speed(2)))
				<< " " << reinmann.rho(0.5*(reinmann.speed(1)+reinmann.speed(2)))
				<< " " << reinmann.rho(0.5*(reinmann.speed(2)+reinmann.speed(3))) << endl;
	}
}

Foam::word exactRiemannSolver::waveToString
(
	const exactRiemannSolver::wave& wave
)
{
	Foam::word waveStr;
	switch(wave)
	{
		case NONE:
			waveStr = "NotAssigned";
		break;
		case SHOCK:
			waveStr = "Shock";
		break;
		case RAREFACTION:
			waveStr = "Rarefaction";
		break;
		case VACUUM:
			waveStr = "Vacuum";
		break;
		case CONTACT:
			waveStr = "ContactDiscontinuity";
		break;
	}
	return waveStr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
