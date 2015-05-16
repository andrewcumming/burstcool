#include <stdio.h>
#include <string.h>
#include "math.h"
#include <stdarg.h>
#include <stdlib.h>
#include <gsl/gsl_sf.h>


#define me 510.999
#define F14 0.352
#define F15 0.648
#define RADa 7.5657e-15
#define PI 3.141592654

#include "../h/odeint.h"
#include "../h/eos.h"
#include "../h/root.h"
#include "../h/vector.h"

void* pt2Object;

// ------------------------ initialise ----------------------------------

void Eos::tidy(void)
{
    free_vector(this->A);
    free_vector(this->Z);
    free_vector(this->X);
}

void Eos::init(int n)
{
  this->ns=n;

  this->A=vector(this->ns);
  this->Z=vector(this->ns);
  this->X=vector(this->ns);
  this->Z2=0.0;
  this->Yn=0.0;
  this->set_Ye=0.0;
  this->set_Yi=0.0;
  this->accr=1;  // default crust composition is HZ with Fe56
  this->gamma_melt=175.0;
  this->Q=900.0; // treat the crust as a liquid for conductivities
  this->B=0.0; // default is unmagnetized 
}

// ------------------------ mean molecular weights ------------------------

double Eos::Yi(void)
  // inverse mean molecular weight per ion
{
  int i;
  double sum=0.0;
  if (set_Yi > 0.0) sum=set_Yi;
  else {
    for (i=1; i<=this->ns; i++)
      sum+=this->X[i]/this->A[i];
  }
  return sum;
}

double Eos::Ye(void)
  // inverse mean molecular weight per electron
{
  int i;
  double sum=0.0;
  if (set_Ye > 0.0) sum=set_Ye;
  else {
    for (i=1; i<=this->ns; i++)
      sum+=this->X[i]*this->Z[i]/this->A[i];
  }
  return sum;
}

double Eos::YZ2(void)
  // sum of (X/A)*Z^2
{
  int i;
  double sum=0.0;
  if (this->Z2 == 0.0) {
    for (i=1; i<=this->ns; i++)
      sum+=this->X[i]*this->Z[i]*this->Z[i]/this->A[i];
  } else {
    sum=this->Z2;  
  }
  return sum;
}

// ------------------------ equation of state ------------------------------

double Eos::pe(void)
  // Calculates the (cgs) electron gas pressure as given by the
  // semi-analytic formula of Paczynski (1983) ApJ 267 315.
{
	double rY, pednr, pedr, pend, ped, ped1,ped2;
 	rY=this->rho*Ye();
	pednr=9.91e-2*pow(rY,5.0/3.0);
	pedr=1.231e1*pow(rY,4.0/3.0);
	ped=1/sqrt((1/pow(pedr,2))+(1/pow(pednr,2)));
  	pend=8.254e-7*1e8*this->T8*rY;
  	return(1e14*sqrt(pow(ped,2)+pow(pend,2)));
}


double Eos::ptot(void)
  // Calculates the total pressure, that is the sum of electron +
  // ion + radiation pressure. Coulomb corrections included if
  // species 1 has Z>2 i.e. not hydrogen or helium, this is a trick to only
  // apply Coulomb corrections in the ocean
{
  double f;
  if (this->Z[1]>2.0) f=1.0+this->Uex()/3.0; else f=1.0;
  double P;
  P=8.254e15*this->rho*this->T8*Yi()*f;   // ions
  P+=pe();                                // electrons  
  P+=RADa*1e32*pow(this->T8,4)/3.0;         // radiation
  return P;
}



double Eos::pemod(void)
  // Modified version of the electron pressure formula
  // to use for heat capacity
  // (Paczynski 1983 ApJ 267 315)
{
  double rY, pednr, pedr, pend, ped;
  rY=this->rho*this->Ye();
  // divide pednr and pedr by appropriate factors
  pednr=9.91e-2*pow(rY,5.0/3.0)/1.32;
  pedr=1.231e1*pow(rY,4.0/3.0)/0.822;
  ped=1/sqrt((1/pow(pedr,2))+(1/pow(pednr,2)));
  pend=8.254e-7*1e8*this->T8*rY;
  return(1e14*sqrt(pow(ped,2)+pow(pend,2)));
}


double Eos::FermiI(int k, double T8, double EF)
  // fitting formula for the generalized Fermi integrals
  // from Chabrier and Potekhin 1998
  // nu=k+1/2, with k=0,1, or 2
{
  double c[3][5]={{0.37045057, 0.41258437, 9.777982e-2, 5.3734153e-3, 3.8746281e-5},
		  {0.39603109, 0.69468795, 0.22322760, 1.5262934e-2, 1.3081939e-4},
		  {0.76934619, 1.7891437, 0.70754974, 5.6755672e-2, 5.5571480e-4}};
  double e[3][5]={{0.43139881, 1.7597537, 4.1044654, 7.7467038, 13.457678},
		  {0.81763176, 2.4723339, 5.1160061, 9.0441465, 15.049882},
		  {1.2558461, 3.2070406, 6.1239082, 10.316126, 16.597079}};
  double x[5]={7.265351e-2, 0.2694608, 0.533122, 0.7868801, 0.9569313};
  double z[5]={0.26356032, 1.4134031, 3.5964258, 7.0858100, 12.640801};
  double h[5]={3.818735e-2, 0.1256732, 0.1986308, 0.1976334, 0.1065420};
  double v[5]={0.29505869, 0.32064856, 7.3915570e-2, 3.6087389e-3, 2.3369894e-5};
  int i;
  double I=0.0, F, R, chi, tau;

  tau=T8/59.4;
  chi=EF/(8.625*T8);

  R=sqrt(chi*(1.0+0.5*chi*tau));

  if (chi*tau > 0.01) {
    F=(chi+1.0/tau)*0.5*R-log(1+tau*chi+sqrt(2.0*tau)*R)/pow(2*tau,1.5);
    if (k>0) F=(2*pow(R,3.0)/3.0-F)/tau;
    if (k>1) F=(2*chi*pow(R,3.0)-5.0*F)/(4.0*tau);
  } else {
    F=pow(chi,1.0*k+1.5)/(1.0*k+1.5);
  }

  if (chi <= 0.6) {
    for (i=0; i<5; i++) {
      I+=c[k][i]*sqrt(1+0.5*e[k][i]*tau)/(exp(-e[k][i])+exp(-chi));
    }
  }
  if (chi > 0.6 && chi < 14.0) {
    for (i=0; i<5; i++) {
      I+=h[i]*pow(x[i],1.0*k)*pow(chi, 1.0*k+1.5)*sqrt(1+0.5*chi*x[i]*tau)/
	  (1+exp(chi*x[i]-chi));
      I+=v[i]*pow(z[i]+chi, 1.0*k+0.5)*sqrt(1+0.5*(z[i]+chi)*tau);
    }
  }
  if (chi >= 14.0) {
    I=F+(PI*PI/6.0)*pow(chi, 1.0*k)*(1.0*k+0.5+0.5*(k+1)*chi*tau)/R;
  }

  return I;
}


double Eos::Fermi_Inv_1_2(double F)
{
  double AN=0.5, RN, DEN, INV, FF;
  int i;
  int M1=2, K1=2, M2=2, K2=2;
  double A1[4]={0.0, 4.4593646e1, 1.1288764e1, 1.0};
  double B1[4]={0.0, 3.9519346e1, -5.7517464, 2.6594291e-1};
  double A2[4]={0.0, 3.4873722e1, -2.6922515e1, 1.0};
  double B2[4]={0.0, 2.6612832e1, -2.0452930e1, 1.1808945e1};
  
  if (F < 4.0) {
    RN=F+A1[M1];
    for (i=M1-1; i>=1; i--) RN=RN*F+A1[i];
    DEN=B1[K1+1];
    for (i=K1; i>=1; i--) DEN=DEN*F+B1[i];
    INV=log(F*RN/DEN);
  } else {
    FF=1.0/pow(F,1.0/(1.0+AN));
    RN=FF+A2[M2];
    for (i=M2-1; i>=1; i--) RN=RN*FF+A2[i];
    DEN=B2[K2+1];
    for (i=K2; i>=1; i--) DEN=DEN*FF+B2[i];
    INV=RN/(DEN*FF);
  }
  
  return INV;
}





double Eos::Chabrier_EF(void)
  // Calculates the Fermi energy in keV including the
  // rest mass using the fit of Chabrier and Potekhin,
  // (1998) Phys Rev E, 58, 4941
  // It uses Antia (1993) to evaluate the inverse Fermi integral
  //
  // rY is (rho/mu_e) in g/cm^3 ;  T is the temperature in K
{
  double EFnr, kT, x, tau, theta;
  double q1,q2,q3, et,etu, corr,F, mc2, rY, T;

  T=this->T8*1e8;
  rY=this->rho*Ye();

  // Electron rest mass in keV
  mc2=510.999;

  // Find kT, x=p_F/m_e c, tau=kT/me and theta=T/T_F
  kT=8.617347*T*1e-8;
  x=1.007e-2*pow(rY,1.0/3.0);
  tau=kT/mc2;
  theta=tau/(sqrt(1.0+x*x)-1.0);

  // Calculate the non-relativistic guess
  F=2.0*pow(theta,-1.5)/3.0;
  EFnr=kT*Fermi_Inv_1_2(F);

  // These functions are defined in CP eq. 24
  if (theta > 69.0) {
    et=1e30; etu=1e-30;
  } else {
    et=exp(theta); etu=1.0/et;
  }
  q1=1.5/(et-1.0);
  q2=12.0+8.0/pow(theta,1.5);
  q3=1.366-(etu+1.612*et)/(6.192*pow(theta,0.0944)*etu+
			   5.535*pow(theta,0.698)*et);

  // This is the correction to the non-relativistic EF
  corr=(1+q1*sqrt(tau)+q2*q3*tau)/(1+q2*tau);
  corr*=tau/(1+(tau/(2*theta)));
  corr=1.5*log(1+corr);

  // return E_F including the rest mass

  // for debugging
  //  printf("Chabrier_EF::::  %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", EFnr, corr, theta, tau,q1,q2,q3,etu,et);

return mc2+EFnr-kT*corr;
}






// ----------------------- thermodynamics --------------------------------

double Eos::chi(double *x)
{
  double x1, x2, p1, p2;
  x1=*x; p1=ptot();
  x2=*x=1.001*x1; p2=ptot();
  *x=x1;
  return (log(p2)-log(p1))/(log(x2)-log(x1));
}


double Eos::CP(void)
  // Calculates specific heat at constant pressure
{
  double cv, chiT, chirho;

  chirho=chi(&this->rho);
  chiT=chi(&this->T8);
  cv=CV();

  return cv+chiT*chiT*ptot()/(this->rho*this->T8*1e8*chirho);
}

double Eos::Gamma1(void)
{
  double chirho, chiT, gam1, cv;

  chirho=chi(&this->rho);
  chiT=chi(&this->T8);
  cv=CV();
  
  gam1=chirho+chiT*chiT*ptot()/(cv*this->rho*1e8*this->T8);

  return gam1;
}  


double Eos::del_ad(void)
  // calculates dlnT/dlnp at constant entropy
{
  double chirho, chiT, gam1, cv;

  chirho=chi(&this->rho);
  chiT=chi(&this->T8);
  cv=CV();
  
  gam1=chirho+chiT*chiT*ptot()/(cv*this->rho*1e8*this->T8);

  return ptot()*chiT/(cv*gam1*this->rho*1e8*this->T8);
}




double Eos::CV(void)
  // Calculates the specific heat at constant volume (density)
{
  double gg, alpha;

  // first the IONS
  gg=this->gamma();

  if (gg < this->gamma_melt) {  // liquid
	    // alpha comes from Chabrier's fit
    double a1,a2,a3,b1,b2,b3,b4;
    a1=-0.9070; a2=0.62954; a3=-0.5*sqrt(3.0)-a1/sqrt(a2);
    b1=4.56e-3; b2=211.6; b3=-1.0e-4; b4=4.62e-3;
    alpha=0.5*pow(gg,1.5)*(a3*(gg-1.0)/pow(gg+1.0,2.0)-a1*a2/pow(gg+a2,1.5))
      +pow(gg,2.0)*(b3*(pow(gg,2.0)-b4)/pow(pow(gg,2.0)+b4,2.0)-
		    b1*b2/pow(gg+b2,2.0));
    cvion=8.3144e7*(1.5+alpha)*this->Yi();
    cv_alpha=alpha;

  } else {  // solid

    double eta=7.76e-5*this->Z[1]*sqrt(this->Yi()*this->rho/(this->A[1]*(1.0-this->Yn)))/this->T8;
 
    double alphaeta=0.399*eta;
    double gameta=0.899*eta;
    double dd,dd1,dd2;
    double x=alphaeta;
    dd1=pow(3.141592654,4.0)/(5.0*pow(x,3.0));
    dd1-=3.0*exp(-x)*(6.0+x*(6.0+x*(3.0+x)))/pow(x,3.0);
    dd2=1.0-0.375*x+0.05*x*x;
    if (dd1 > dd2) dd=dd2; else dd=dd1;
    cvion=8.3144e7*this->Yi()*(8.0*dd-6*alphaeta/(exp(alphaeta)-1.0)+(pow(gameta,2.0)*exp(gameta)/pow(exp(gameta)-1.0,2.0)));
    
  }

  // RADIATION
   	//cvrad=3.0256e10*pow(this->T8,3)/this->rho;
    cvrad=4.0*RADa*1e24*pow(this->T8,3)/this->rho;
//  cvrad=0.0;
  
  // ELECTRONS
  { // modified version of Paczynksi's fit for cve
    double dT,temp,p1,p2;
    temp=1.01*this->T8; dT=temp-this->T8;
  // p1=this->pe(); this->T8+=dT; p2=this->pe(); this->T8-=dT;
      p1=this->pemod(); this->T8+=dT; p2=this->pemod(); this->T8-=dT;
    cve=(1/((this->f()-1)*this->rho))*1e-8*(p2-p1)/dT;
  }

 

  // total
return cve+cvion+cvrad;
}

/*
double Eos::CV(void)
  // Calculates the specific heat at constant volume (density) by 
  // performing a numerical differentiation. The electron specific heat
  // is calculated using the fitting formula given by Paczynski (1983).
{
  double cv, T1, T2, p1, p2, U1, U2;
  T1=this->T8; p1=pe(); U1=this->Uex();
  this->T8=T2=1.001*T1; p2=pe(); U2=this->Uex();
  this->T8=T1;
  cv=0.0;
  // ions
  cv+=8.3144e7*1.5*this->Yi();
  if (this->Z[1] > 1.0) 
    cv+=8.3144e7*this->Yi()*(U2*T2-U1*T1)/(T2-T1);
  // electrons
  cv+=(1/((f()-1)*this->rho))*1e-8*(p2-p1)/(T2-T1);
  // radiation
  cv+=3.0256e10*pow(T1,3)/this->rho;
  return cv;
}
*/

double Eos::f(void)
  // Calculates f = dln ped/dln rho, using the fitting formula given
  // by Paczynski (1983).
{
  double rY, pednr, pedr, ped;
  rY=this->rho*Ye();
  pednr=9.91e-2*pow(rY,5.0/3.0);
  pedr=1.231e1*pow(rY,4.0/3.0);
  ped=1/sqrt((1/pow(pedr,2))+(1/pow(pednr,2)));
  return (5.0*pow(ped/pednr,2) + 4.0*pow(ped/pedr,2))/3.0;
}



// --------------------------- opacity ------------------------------------


double Eos::eps_nu(void)
  // Calculates neutrino emissivity (erg/g/s)
  // by plasma process from Schinder et al. 1987
{
	double a0, a1, a2, b1, b2, b3, c;
	double xi, xi2, xi3, la, la2, la3, g, K;

	// variables
	la=this->T8/59.302; la2=la*la; la3=la2*la;
	xi=pow(this->rho*this->Ye()*1e-9, 1.0/3.0)/la;
	xi2=xi*xi; xi3=xi2*xi;

	// 1. plasma

	// these coefficients valid for 10^8<T<10^11 K
	a0=2.146e-7; a1=7.814e-8; a2=1.653e-8;
	b1=2.581e-2; b2=1.734e-2; b3=6.990e-4;
	c=0.56457;
	K=pow(this->rho*this->Ye(),3.0);
	// formula from Schinder et al.
	Q1=K*exp(-c*xi)*(a0+a1*xi+a2*xi2)/(xi3+(b1/la)+(b2/la2)+(b3/la3));

	// 2. pair
	// coefficients valid for 10^8 < T < 10^11 K
	a0=5.026e19; a1=1.745e20; a2=1.568e21;
	if (this->T8 < 100.0) {  // 10^8<T<10^10 K
		b1=9.383e-1; b2=-4.141e-1; b3=5.829e-2;
		c=5.5924;
	} else { // 10^10 < T < 10^11 K
		b1=1.2383; b2=-8.141e-1; b3=0.0;
		c=4.9924;
	}
	g=1.0-13.04*la2+133.5*la2*la2+1534*la2*la2*la2+918.6*la2*la2*la2*la2;
	K=g*exp(-2.0/la);
	double qpair;
	qpair=pow(10.7480*la2+0.3967*sqrt(la)+1.0050,-1.0)
		* pow(1.0 + this->rho*this->Ye()/(7.692e7*la3+9.715e6*sqrt(la)),-0.3);
	// formula from Schinder et al.
	Q2=(1.0+0.10437*qpair)*K*exp(-c*xi)*(a0+a1*xi+a2*xi2)/(xi3+(b1/la)+(b2/la2)+(b3/la3));

	// 3. Brems    formula from Haensel et al. 96
	//Q3=0.3229*this->rho*this->YZ2()*pow(this->T8,6.0);
	// should now multiply by a order unity factor, add this later

	if (this->gamma() < this->gamma_melt && this->rho <1e10) {
		double L;
		double A, B, eta, t, Z;
		Z=this->Ye()/this->Yi();
	// t = kT/2p_Fc (see their eq. 24)
		t=this->T8/(118.6*(sqrt(1.0+pow(this->x(),2.0))-1.0));

	// finite core radius : value depends on whether above or
	// below neutron drip
		if (this->Yn == 0.0) eta=0.16*pow(this->rho*1e-12,1.0/3.0);
		else eta=0.25*pow(this->rho*1e-12*this->Ye(),1.0/3.0);
	//    if (eta > 0.2) eta = 0.2;

		A=0.269+20.0*t+0.0168*Z+0.00121*eta-0.0356*Z*eta+0.0137*Z*Z*t+1.54*Z*t*eta;
		B=1.0+180.0*t*t+0.483*t*Z+20.0*t*Z*eta*eta+4.31e-5*Z*Z;

		L=A/pow(B,0.75);
	//L=1.0;

		Q3=L*0.3229*this->rho*this->YZ2()*pow(this->T8,6.0);

	//Q3=L*0.3229*this->rho*pow(this->T8,6.0);
	//Q3*=1.0-this->Yn;  // mass fraction of nuclei
	//Q3*=this->Z[1]*this->Z[1]/this->A[1];
	//   Q3=0.0;

	} else {

	// Fit from Kaminker et al. 99
		double r=log10(this->rho*1e-12);
		double t=log10(this->T8);
		double r0=2.8e14;

		Q3=11.204 + 7.304*t + 0.2976*r - 0.370*t*t + 0.188*t*r - 0.103*r*r
			+ 0.0547*t*t*r - 6.77*log10(1+0.228*this->rho/r0);

		if (this->accr) {
	// Dany suggested the following adjustments to approximate accreted matter
			if (r < 0.1) Q3-=0.2;
			else {
				if (r < 1.0) Q3-=0.3;
				else Q3-=0.4;
			}
		}

		Q3=pow(10.0,Q3);
	}


//	Q3=0.0;
	return (Q1+Q2+Q3)/this->rho;
}





double Eos::opac(void)
  // Calculates the opacity
{
	double kappa, ef, eta, TP;
	int i;

	// Fitting formula for the electron scattering opacity from Paczynski 1983 ApJ
	this->kes=(0.4*Ye())/((1+2.7e11*this->rho*pow(1e8*this->T8,-2.0))*
		(1+pow(this->T8/4.5,0.86)));

	// Fermi energy
	ef=Chabrier_EF();
	if (ef == 0) {
		eta=Fermi_Inv_1_2(1.105e-4*this->rho*Ye()/pow(this->T8,1.5));
		ef=eta*8.617*this->T8;
	} else {
		eta=(ef-me)/(8.617*this->T8);
	}

	// Free-Free opacity
	this->kff=7.53e-6*this->rho*Ye()/pow(this->T8, 3.5);

	kgaunt=0.0;
	for (i=1; i<=this->ns; i++)
		kgaunt+=this->Z[i]*this->Z[i]*this->X[i]*gff(this->Z[i],eta)/this->A[i];
	this->kff*=kgaunt;

	// total radiative opacity
	this->kappa_rad=this->kff+this->kes;

	// "non-additivity" factor from Potekhin et al. (2001) eqs 19-20
	double f=this->kff/this->kappa_rad;
	double TRy = 100.0*this->T8/(0.15789*this->YZ2()/this->Yi());
	double A = 1.0 + (1.097+0.777*TRy)*pow(f,0.617)*pow(1.0-f,0.77)/(1.0+0.536*TRy);
	// or use my version:
	//double A = 1.0+0.42*exp(-pow(0.37*log(this->kff/(0.8*this->kes)),2.0));
	this->kappa_rad*=A;

	// add correction for frequencies < plasma frequency
	// (not sure where this came from ??) 
	//TP=3.3e-3*sqrt(this->rho*this->Ye());  // plasma temperature
	//kappa_rad*=exp(TP/this->T8);
	//double up=TP/this->T8;
	//kappa_rad*=(1-2.448*pow(0.1*up,2.0)+16.40*pow(0.1*up,6.0));

	// Conduction
	this->kcond=3.024e20*pow(this->T8,3)/(K_cond(ef)*this->rho);

	// Add up opacities  
	kappa=1.0/((1.0/this->kcond)+(1.0/this->kappa_rad));

	return kappa;
}


double Eos::gff(double Z1, double eta)
	// Calculates the free-free Gaunt factor for element with charge Z1 
	// using a fitting formula described in appendix of Schatz et al (1999) ApJ
	// Agrees to 10% with Itoh et al. 1991
{
  double gaunt, x, rY, T8_32, gam, ef;

  rY=this->rho*Ye();
  T8_32=pow(this->T8, 1.5);

  if (eta < 100.0) x=log(1.0+exp(eta)); 
  else x=eta; // make sure it doesn't freak out for extremely large eta

  // normalisation and degeneracy piece
  gaunt=1.16*8.02e3*x*T8_32/rY;

  x=pow(1+x,2.0/3.0);
  gam=sqrt(1.58e-3/this->T8)*Z1;

  // Elwert factor
  gaunt*=(1.0-exp(-2*PI*gam/sqrt(x+10.0)))/(1.0-exp(-2*PI*gam/sqrt(x)));

  // relativistic piece
  gaunt*=1.0+pow(this->T8/7.7, 1.5);

  // send it back
  return gaunt;
}


double Eos::Uex(void)
  // Coulomb correction Uex/kT
{
  double u,g,g2,g3,g14;
  g=this->gamma(); g2=g*g; g3=g2*g; g14=pow(g,0.25);

  if (g < this->gamma_melt) {
    u=-0.89813*g+0.98686*g14-0.91095+0.25098/g14;
  } else {
    u=-0.89593*g+1.5+9.65/g+840/g2+1.101e5/g3;
  }

  return u;
}


double Eos::Fep(int flag)
  // "Coulomb log" for electron-phonon scattering 
  // (Baiko & Yakovlev 1995,1996)
  // if flag=0 electrical conductivity; flag>0 thermal conductivity
{
  double R0, R1, R2, G0, G2, t, u2, u1, s, F, K0,K1;
  double alpha, alpha0, a0, a2, x, beta, AA, ZZ;
  double P0,g, K2,P2, c1, c2;

  // constants -- use values for bcc crystal
  a0=0.0174; a2=0.0118; u2=13.0; u1=2.8;

  AA=this->A[1]; ZZ=this->Z[1];
  x=this->x(); beta=x/sqrt(1+x*x);
  t=0.804*this->T8*(0.5*AA/ZZ)/sqrt(1e-9*this->rho);
  s=pow(4*ZZ,-2.0/3.0)+2.323e-3/beta;
  alpha0=1.683*sqrt(x/(AA*ZZ));
  alpha=alpha0*(0.5*u1*exp(-9.1*t)+t*u2);
  //alpha=1e-6;  //  small alpha is the Yakovlev & Urpin limit
  
  G0=u2*t/sqrt(t*t+a0); 
  R0=(exp(-alpha*s)-exp(-alpha))/alpha;
  R1=2*(exp(-alpha*s)*(1+alpha*s)-exp(-alpha)*(1+alpha))/
    (alpha*alpha);
  K0=2*R0-beta*beta*R1;
  
  // correction for finite nuclear size
  if (this->rho < 4e11) g=0.16*pow(this->rho*1e-12,1.0/3.0);
  else g=0.25*pow(this->rho*1e-12*this->Ye(),1.0/3.0);
  //g=0.0;  //switch off finite size effects
  P0=4.787-0.0346*ZZ;
  R2=(exp(-alpha*s)*(alpha*alpha*s*s+2*alpha*s+2)-exp(-alpha)*(alpha*alpha+2*alpha+2))/(alpha*alpha*alpha);
  c1=pow(1.0+pow(18.0*ZZ*PI,2.0/3.0)*g*g*(0.5*R1-beta*beta*R2)/(2.5*K0*P0),-P0);
  
  F=G0*K0*c1;
  
  if (flag > 0) { // thermal conductivity so add an extra piece
    P2=2.729-0.0204*ZZ;
    R2=this->Eep(alpha*s)-this->Eep(alpha);
    G2=t/(PI*PI*pow(t*t+a2,1.5));
    K2=0.5*R2-0.5*beta*beta*R0;
    // correction for finite nuclear size
    c2=pow(1.0+pow(18.0*PI*ZZ,2.0/3.0)*g*g*0.5*K0/(10.0*K2*P2),-P2);
    F+=G2*(3*K2-0.5*K0)*c2;
  }
  return F;
}

double Eos::Eep(double q)
  // used by Fep() to calculated thermal conductivity piece
  // Baiko & Yakovlev 1995
{
  double q2,q3,q4,qu;
  q2=q*q; q3=q2*q; q4=q3*q; qu=1.0/q;
  return exp(-q4/(q3+0.1397))*(log(1+qu)-0.5772/(1+2.2757*q2));
}


double Eos::lamei(int n)
// n=1 for thermal n=0 for electrical
{

  if (this->T8 < 0.0 || isnan(this->T8)) return 1.0;

//if (this->rho > 3e14) return 1.0;

  if (0) {
  /* Yakovlev & Urpin */

  double x,x1,x2,lam;
  x1=this->x();
  x2=0.22*sqrt(this->T8);
  if (x1>x2) x=x1; else x=x2;
  lam=127*x*sqrt((3.0/this->gamma())+1.5)/
    pow(this->rho*this->Yi(),1.0/3.0);
  return log(lam)-0.5*x*x/(1+x*x);
  
  } else {

  // Potekhin et al. 1999

  double G, eta, eta0, beta, x, vc, gam;
  double ZZ=this->Ye()/this->Yi();
  
	x = this->x();
	eta=this->T8/(0.07832*sqrt(this->rho*1e-6)*this->Ye());
  	eta*=sqrt(1.0-this->Yn);
 
//	if (this->rho < 7.09e3*pow(1e-12*this->B,1.5)/this->Ye()) {
//		double xr = 2.96e-5*this->rho*this->Ye()*1e12/this->B;			
//		if (xr > x) x=xr;
//	}
//			if (xr < sqrt(1.38e-8*this->T8/(9.11e-28*9e20))) xr=1e-10;
//		eta = 1.38d-8*this->T8/(9.11e-28*9e20*(sqrt(x*x+1.0)-1.0));
//	}
vc=x/sqrt(1+x*x);

  gam=this->gamma();

//  printf("T8=%lg, rho=%lg, x=%lg gam=%lg\n", this->T8, this->rho,x,gam);

  //if (gam < 1.0) return 1.0;

  beta=ZZ*PI*vc/137.0;
  eta0=0.19/pow(ZZ,1.0/6.0);
  
  // the easy part is to calculate G
  G=eta*(1.0+0.122*beta*beta)/sqrt(eta*eta+eta0*eta0);
  if (n==1) G+=0.0105*(1.0-1.0/ZZ)*(1+beta*vc*vc*vc)*eta/pow(eta*eta+0.0081,1.5);
  
  // next the s and w parameters
  double rtf, rd;
  rtf=1.0/(PI*137.0*beta);
  rd=(59.41/this->T8)*ZZ*x/(PI*3.0*137.0);

  // printf("%lg %lg %lg %lg %lg %lg ", this->T8, this->rho, beta, rd, gam, x);
  s=exp(-beta)*(rtf+rd*(1.0+0.06*gam)*exp(-sqrt(gam)));
  w=13.0*(1.0+beta/3.0)/rd;

  if (isnan(s) || isnan(w)) printf("%lg %lg %lg\n", this->T8,s, w);
    
  // now calculate the exponential integrals    
  //  double L1, L2;

 // if (1) {

 	if (w > 50.0) {

    L1=0.5*(log((1.0+s)/s)-1.0/(1.0+s));
    L2=(2.0*s+1.0)/(2.0*s+2.0)-s*log((1.0+s)/s);

  } else {

  if (s < 1e-2 && s < 1e-2/w) {
    L1=0.5*(expint(1,w)+log(w)+0.5772);
    L2=(exp(-w)-1.0+w)/(2.0*w);
    
  } else {
    double I1,I2;
      // printf("s=%lg w=%lg    rho=%lg T8=%lg\n",s,w,this->rho, this->T8);
    I1=expint(1,s*w);
    //printf("I1=%lg\n",I1);
    I2=expint(1,w*(1.0+s));
    //printf("I2=%lg\n",I2);
 
    L1=log((1.0+s)/s)+(s/(1.0+s))*(1.0-exp(-w))-(1.0+s*w)*exp(s*w)*
      (I1-I2);
    L1/=2.0;
    
    L2=((exp(-w)-1.0+w)/w)-(s*s/(1.0+s))*(1.0-exp(-w))-2*s*log((1.0+s)/s)+
      s*(2.0+s*w)*exp(s*w)*(I1-I2);
    L2/=2.0;
    
    // printf("L1=%lg L2=%lg I1-I2=%lg\n", L1,L2,I1-I2);
  }
  }
  // put it all together
  G*=(L1-vc*vc*L2);

  double D=exp(-0.42*sqrt(x/((this->A[1]*(1.0-this->Yn))*this->Z[1]))*3.0*exp(-9.1*eta));
  G*=D;

  return G;

  }  

}  

double Eos::expint(int n, double x)
{
	if (n != 1) printf("I only know how to do Ei1(x)!\n");
	return gsl_sf_expint_E1(x);
}


double Eos::find_rho(void)
{
  pt2Object=(void*) this;
if (Wrapper_find_rho_eqn(1e-6) > 0.0) return 1e-6;
  return zbrent(Wrapper_find_rho_eqn,1e-6,1e15,1e-6);
}


double Eos::Wrapper_find_rho_eqn(double r)
{
  Eos* mySelf = (Eos*) pt2Object;
  return mySelf->find_rho_eqn(r);
}

double Eos::find_rho_eqn(double r)
{
  this->rho=r;
  return this->ptot()-this->P;
}

double Eos::x(void)
{
  double x; 
  x=pow(this->Chabrier_EF()/511.0,2)-1.0; if (x<0.0) x=1e-10;   
  x=sqrt(x);
  return x;
}

double Eos::eta(void)
{
  return (this->Chabrier_EF()-511.0)/(8.625*this->T8);
}

double Eos::gamma(void)
{ 
//  return 0.11*(this->YZ2()/this->Yi())*pow(this->rho*1e-5*Yi(),1.0/3.0)/this->T8;
  return 0.11*this->Z[1]*this->Z[1]*pow(this->rho*1e-5*Yi(),1.0/3.0)/this->T8;
}

double Eos::K_cond(double ef)
  // Calculates the conductivity due to electron-ion and
  // electron-electron collisions
  // ef is the Fermi energy in keV
{
  double x, lam, y, y3, rY, x2, K;
  double gam, f_c, theta, beta, corr;

  // set up parameters
  rY=this->rho*this->Ye();
  x=this->x();

  //double xx1=this->x();
  //double xx2=0.26*sqrt(this->T8);
  // if (xx1>xx2) x=xx1; else x=xx2;
  
  x2=sqrt(1+x*x); beta=x/x2;
  gam=this->gamma();

  // This is the Coulomb logarithm from Yakovlev and Urpin
 // lam=log(pow(2*PI*Ye()/(3*Yi()),1.0/3.0)*sqrt(1.5+3.0/gam));
 // lam-=0.5*beta*beta;
//  this->lambda2=lam;

  // get the Coulomb logarithm for electron-ion collisions
  // from Potekhin et al. 1999
	if (this->rho > 1e3)  this->lambda2=this->lamei(1); else this->lambda2=1.0;

  // electron-electron collisions
  // Note that Potekhin et al 1997 which is where we get the J(x,y)
  // function has a misprint in the prefactor for f_ee. 
  // The correct expression is in Timmes 1992 or in 
  // Potekhin et al. (1999).
  
  y=5.771e-3*sqrt(rY/x2)/this->T8;
  f_ee=5.11e15*this->T8*this->T8*pow(x,1.5)*J(x,y)/pow(1+x*x,1.25);
  
  if (gam < this->gamma_melt || this->Q == 900.0) { // if Q=900 treat as liquid 

    // The electron-ion collision frequency
    f_ei=1.76e16*this->lambda2*x2*YZ2()/Ye();

    // The collision frequencies add
    f_c=f_ee+f_ei;

  } else { // solid --- NB assumes A=2Z and single species

    // The electron-ion collision frequency (ie. phonons)
    // is calculated as given by Potekhin et al. 1999
    f_ep=1.76e16*this->lambda2*x2*YZ2()/Ye();
    
    // add exponential suppression when the Umklapp scatterings freeze out
    {
      double TU=2.2e8*sqrt(1e-12*this->rho)*this->Ye()*pow(this->Z[1]/60.0,1.0/3.0);
      if (this->T8<1e-8*TU) 
	f_ep*=exp(-1e-8*TU/this->T8);
    }

    //   f_ep*=(1.0-this->Yn);  // to agree with Ed

    /* old phonons from Yakovlev & Urpin
       theta=0.56*sqrt(1e-9*this->rho)/this->T8;
       lam=(2-beta*beta)/(beta*sqrt(1+pow(theta/3.5,2.0)));
       lam+=pow(theta/5.1,2.0)*(3*this->lambda2-1+0.5*beta*beta)/
       (beta*pow(1+pow(theta/4.2,2.0),1.5));
       f_c=1.24e18*this->T8*lam;
    */
    /* and from Baiko & Yakovlev 
       f_ep=9.55e16*this->T8*this->Fep(1)/beta; // phonons
    */

    // Impurity scattering
    // Coulomb log from Itoh & Kohyama 1993
    {
      double ka, sm1;
		// the following is eq.(20) of IK93 for ka
      //ka=1.92*pow(this->Ye()/this->Yi(),1.0/3.0);
	  // instead, we use the substitution  ka -> 2k/kTF and kTF is given by eq.(3) of
    	// Potekhin et al. 1999
      ka=sqrt(137.0*PI*beta);

	  // and then put into eqs (10,16,17) of IK93
      sm1=0.5*log(1.0+0.4*ka*ka);
      lam=sm1*(1.0+2.5*beta*beta/(ka*ka))-0.5*beta*beta;
    }      

    if (this->Yn > 0.0) f_eQ=1.76e16*this->Q*lam*x2/this->Z[1];
    else f_eQ=1.76e16*this->Q*lam*x2/this->Z[1];
 
    /*
    // multiply by two to get Ed's answer
    f_eQ*=2;
    if (this->rho > 1e13) f_eQ*=1.4;  
    // and by a further 40% in the inner crust
    else {
      if (this->rho > 2e12) f_eQ*=1.2; 
    }
    */
    // sum of phonons and impurities and electrons
    f_c=f_eQ+f_ep+f_ee;
  }
  
  // the conductivity is then as given by Yakovlev & Urpin
  K = 4.116e27*this->T8*rY/(x2*f_c);

  // correction due to thermoelectric field
  //corr=9.34e-4*pow(this->T8,2.0)*(1+x*x)/pow(x,4.0);
  //corr*=pow(6.0-2.0*beta*beta-(1.0-beta*beta+beta*beta*beta*beta)/lam,2.0);
  //K=K*(1.0-corr); if (K<0.0) return 1e-10;

  return K; 
  
}



extern "C"{
  void condegin_(double *temp,double *densi,double *B,double *Zion,double *CMI,double *CMI1,double *Zimp, double *RSIGMA,double *RTSIGMA,double *RHSIGMA,double *RKAPPA,double *RTKAPPA,double *RHKAPPA);
   }

double Eos::potek_cond()
{
	// This is a wrapper for Potekhin's conductivity routine
	// and returns the non-magnetic thermal conductivity in cgs units
	double s1,s2,s3,k1,k2,k3;
	double null=0.0, Zimp=sqrt(this->Q), AA=this->A[1]*(1.0-this->Yn);
	double temp=this->T8*1e2/5930.0;
	double rr=this->rho/(this->A[1]*15819.4*1822.9);
	condegin_(&temp,&rr,&null,&this->Z[1],&AA,&this->A[1],&Zimp, &s1,&s2,&s3,&k1,&k2,&k3);
	return k1*2.778e15;
}



double Eos::J(double x,double y)
{
  // from Potekhin, Chabrier, & Yakovlev 1997
  double x2=x*x;
  double b2=x2/(1.+x2);
  double y3=y*y*y;
  double y4=y3*y;
  return (1.+0.4*(3.+1./x2)/x2)*(y3*pow(1.+0.07414*y,-3.)*
                                 log((2.810-0.810*b2+y)/y)/3.+
                                 pow(PI,5)*y4*pow(13.91+y,-4.)/6);
}

