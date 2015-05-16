// odeint.cc
//

#include "../h/vector.h"
#include <stdio.h>
#include "math.h"
#include "../h/odeint.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>


void *pointer_to_OdeInt;

int Ode_Int::gsl_derivs (double t, const double y[], double dydt[], void * params)
{
	Ode_Int *myself = (Ode_Int*) pointer_to_OdeInt;

	for (int i=1; i<=myself->nvar; i++)
		myself->derivs_y[i] = y[i-1];

	myself->delegate->derivs(t,myself->derivs_y,myself->derivs_dydt);

	for (int i=1; i<=myself->nvar; i++)
		dydt[i-1] = myself->derivs_dydt[i];	

	return GSL_SUCCESS;
}


int Ode_Int::gsl_jacobn (double t, const double y[], double *dfdy,double dfdt[], void *params)
{
	Ode_Int *myself = (Ode_Int*) pointer_to_OdeInt;

	for (int i=1; i<=myself->nvar; i++)
		myself->derivs_y[i] = y[i-1];
		
	myself->delegate->derivs(t,myself->derivs_y,myself->derivs_dydt);
	
	myself->delegate->jacobn(t,myself->derivs_y,myself->derivs_dydt,myself->derivs_dfdy,myself->nvar);

	for (int i=1; i<=myself->nvar; i++)
		for (int j=1; j<=myself->nvar; j++)
			dfdy[(i-1)*myself->nvar + (j-1)] = myself->derivs_dfdy[j][i];

	return GSL_SUCCESS;	
}

void Ode_Int::tidy(void)
{
	free_vector(this->ystart);
	free_vector(this->ynext);
	free_vector(this->hstr);
	free_vector(this->xp);
	free_matrix(this->yp,this->nvar,this->kmax);
	free_matrix(this->dydxp,this->nvar,this->kmax);
}

void Ode_Int::init(int n, Ode_Int_Delegate *delegate)
{
	this->delegate = delegate;

	this->kmax=900000;
	this->nvar=n;
	this->ignore=0;
	this->dxsav=0.0;
	this->minstep=0.0;
	this->hmax=1e12;

	this->xp=vector(this->kmax);
	this->yp=matrix(this->nvar,this->kmax);
	this->dydxp=matrix(this->nvar,this->kmax);
	this->hstr=vector(this->kmax);
	this->ystart=vector(this->nvar);
	this->ynext=vector(this->nvar);

	this->stiff = 0; // default is non-stiff eqns.
	this->verbose = 0; // turn off output
	this->tri=0; // don't assume a tridiagonal Jacobian
		// if we are using GSL, this has no effect

	this->use_gsl=1;  // use the GSL integrator
}

void Ode_Int::set_bc(int n, double num)
{
  this->ystart[n]=num;
}

double Ode_Int::get_d(int n, int i)
{
  return this->dydxp[n][i];
}

double Ode_Int::get_x(int i)
{
  return this->xp[i];
}

double Ode_Int::get_y(int n, int i)
{
  return this->yp[n][i];
}

void Ode_Int::go(double x1, double x2, double xstep, double eps)
{
	go_gsl(x1,x2,(long int)((x2-x1)/xstep),eps,1);
}



void Ode_Int::go_simple(double x1, double x2, int nstep)
{
	go_gsl(x1,x2,nstep,1e-6,0);
}


void Ode_Int::go_gsl(double x1, double x2, long int nsteps, double eps, int log_flag)
{
	pointer_to_OdeInt = (void*) this;

	if (this->verbose) printf("Number of steps=%ld\n",nsteps);

	double xstep;
	if (log_flag) xstep = pow(10.0,log10(x2)/(1.0*nsteps));
	else xstep = (x2-x1)/(1.0*nsteps);

	this->derivs_y=vector(this->nvar);
	this->derivs_dydt=vector(this->nvar);
	this->derivs_dfdy=matrix(this->nvar,this->nvar);
	for (int i=1; i<=this->nvar; i++)
		for (int j=1; j<=this->nvar; j++)
			this->derivs_dfdy[i][j] = 0.0;

	this->sys = (gsl_odeiv2_system) {gsl_derivs,gsl_jacobn,this->nvar,NULL};
//	if (this->stiff) this->step=gsl_odeiv2_step_alloc (gsl_odeiv2_step_bsimp,this->nvar);
//	else this->step=gsl_odeiv2_step_alloc (gsl_odeiv2_step_rkf45,this->nvar);
//	this->control=gsl_odeiv2_control_y_new(0.0,eps);
//	this->evolve=gsl_odeiv2_evolve_alloc(this->nvar);
	this->driver=gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msbdf,xstep,0.0,eps);
	
	double x = x1;
	double h = xstep;
	
	for (int i=1; i<=this->nvar; i++) {
		this->ynext[i-1]=this->ystart[i];
		this->yp[i][1]=this->ystart[i];
	}	
	this->xp[1]=x1;
	this->kount=1;
	
	int status;

	for (int j=1; j<=nsteps; j++) {
	
		double xnext;
		if (log_flag) xnext = pow(10.0,log10(x2)*j/(1.0*nsteps));
		else xnext = (x2-x1)*j/(1.0*nsteps);
		
		status=gsl_odeiv2_driver_apply (this->driver,&x,xnext,this->ynext);

		if (this->verbose) printf("%lg %lg %lg %d\n",x1,x2,xnext,status);

		if (status != GSL_SUCCESS) break;

		this->kount++;
		if (this->kount == this->kmax) {
			printf("Maximum number of steps reached! Stopping integrator.\n");
			break;
		}

		this->xp[this->kount]=x;
		for (int i=1; i<=this->nvar; i++) {
			this->yp[i][this->kount]=this->ynext[i-1];
			this->dydxp[i][this->kount]=0.0;
		}
		
	}


/*
	while (x < x2) {
	//	status = gsl_odeiv2_evolve_apply(this->evolve, this->control, this->step,
	//			&this->sys, &x, x2, &h, this->ynext);
			status=gsl_odeiv2_driver_apply (this->driver,&x,x2,this->ynext);
		
		if (status != GSL_SUCCESS) break;
		this->kount++;
		if (this->kount == this->kmax) {
			printf("Maximum number of steps reached! Stopping integrator.\n");
			break;
		}
		this->xp[this->kount]=x;
		for (int i=1; i<=this->nvar; i++) {
			this->yp[i][this->kount]=this->ynext[i-1];
			this->dydxp[i][this->kount]=0.0;
		}
	}
	*/

//	gsl_odeiv2_evolve_free(this->evolve);
//	gsl_odeiv2_control_free(this->control);
//	gsl_odeiv2_step_free(this->step);
	gsl_odeiv2_driver_free (this->driver);
	
	free_matrix(this->derivs_dfdy,this->nvar,this->nvar);
	free_vector(this->derivs_y);
	free_vector(this->derivs_dydt);
}
