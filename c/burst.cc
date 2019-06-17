#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "../h/vector.h"
#include "../h/root.h"
#include "../h/ns.h"
#include "../h/timer.h"
#include "../h/burst.h"

void *pt2Burst;

Burst::Burst() {
	// initialize default parameters
	this->mass = 1.4;
	this->radius = 12.0;
  	this->temperature_slope=-1;
	this->E18=0.2;
	this->yb=1e12;
	this->ydeep_factor=100.0;
	this->yt=1e8;
	this->convection_flag=1;
	this->time_to_run=1e5;

	this->deep_composition_flag=1;
	this->shallow_composition_flag=1;
	
	// to model URCA cooling:
	this->icool = 32;  // index of the grid cell for the cooling source
	this->L34 = 0.0;   // luminosity of the cooling source

	// Can use this to turn off neutrino emission
	this->nuflag=1;
	// If outer_boundary_flag is set, use the outer grid from makegrid, otherwise
	// if the flag ==0 assume F\propto T^4 on the outer boundary
	this->outer_boundary_flag=1;

	// Output to files?
	this->output = 1;

	// Envelope files
	this->envelope_file = 'none';
	this->env_g = 2.45e14;

}


void Burst::setup(void) {

	if (this->yb < 1000.0) this->yb=pow(10.0,this->yb);
	if (this->yt < 1000.0) this->yt=pow(10.0,this->yt);

	set_ns_parameters(this->mass,this->radius,&this->g,&this->ZZ);

  	// initialize the numerical grid
 	// the last parameter is composition
 	// 1=Ni56  2=Ca40  3=28Si  4=12C
  	set_up_grid(50);

	// read in Tb-this->TEFF relation
	// first param = He column;  second 0=Fe, 1=Si
	if (this->outer_boundary_flag) get_TbTeff_relation(4,0);

	// precalculate density, CP, K etc as a function of T on the grid
	clock_t timer;
	start_timing(&timer);
	precalculate_vars();
	stop_timing(&timer,"precalculate_vars");

	// initialize the integrator
	this->ODE.init(this->N+1,dynamic_cast<Ode_Int_Delegate *>(this));
	this->ODE.verbose=0;   // print out each timestep if this is set to 1
	this->ODE.stiff=1; this->ODE.tri=1;  // stiff integrator with tridiagonal solver

	// initial temperature profile
	pt2Burst = (void *)this;

	if (this->convection_flag) {   // adiabatic or power law temperature profile
		double beta = 0.99;
		double T2=pow(3.0*beta*this->yb*this->g/7.5657e-15,0.25);	
		double Tb=zbrent(Wrapper_set_initial_temperature_profile,5e8,T2,1e-3);
		printf("initial Tb=%lg \n",Tb);
	} else {    // instantaneous burn
		for (int i=this->N+1; i>=1; i--) {
			double Tf, Ti=2e8;
			if (this->y[i] <= this->yb+1.0) {	
				this->heat_y = this->y[i];
				this->heat_Ti = 2e8;
				double beta = 0.5;
				double T2=pow(3.0*beta*this->y[i]*this->g/7.5657e-15,0.25);	
				Tf = zbrent(Wrapper_heat,2e8,T2,1e-3);
			} else {
				Tf = Ti;
			}
			this->ODE.set_bc(i,Tf);
		}
	}

	// open output files
	if (this->output) {
		this->fp=fopen("out/out","w");
		this->fp2=fopen("out/prof","w");
		this->fp3=fopen("out/summary","a");
		fprintf(this->fp,"%d\n",this->N);
	}
}

	
Burst::~Burst() {
	// destructor
  	if (this->output) {
		fclose(this->fp);
  		fclose(this->fp2);
		fclose(this->fp3);
	}
  	this->ODE.tidy(); 
  	free_vector(this->CP);
  	free_vector(this->y);
  	free_vector(this->K);
  	free_vector(this->F);
  	free_vector(this->NU);
}

double Burst::Wrapper_set_initial_temperature_profile(double Tb)
{
	Burst *myself = (Burst *) pt2Burst;
	return myself->set_initial_temperature_profile(Tb);
}


double Burst::Wrapper_heat(double Tf)
{
	Burst *myself = (Burst *) pt2Burst;
	return myself->heat(Tf) - myself->E18*1e18;
}



double Burst::set_initial_temperature_profile(double Tb)
{
	// initial temperature profile
	double Ti = 2e8, Tf=Tb, energy = 0.0;
	for (int i=this->N+1; i>=1; i--) {
		if (this->y[i] <= this->yb+1.0) {
	//			this->ODE.set_bc(i,Tf*pow(this->y[i]/this->yb,this->temperature_slope)); 
			double T2 = Tf;
			double Prad=2.521967e-15*pow(Tf,4);
  			if (Prad > this->y[i]*this->g) T2*=pow(this->y[i]*this->g/Prad,0.25)*0.99;
			this->ODE.set_bc(i,T2);
			this->heat_y = this->y[i];
			this->heat_Ti = Ti;
			energy += heat(T2)*this->y[i]*this->dx;
			if (this->temperature_slope < 0.0) {  // isentropic temperature profile
				set_composition(this->shallow_composition_flag);
				this->EOS.T8 = 1e-8*T2;
				this->EOS.rho = this->EOS.find_rho();
				double del_ad = this->EOS.del_ad();
				Tf *= (1.0 - this->dx*del_ad); 
	//			printf("%lg %lg %lg %lg\n",this->y[i],Tf,T2,del_ad);						
				
			} else {    // power law temperature profile
			//printf("%lg %lg %lg %lg\n",this->y[i],Tf,T2,heat(Ti,T2,this->y[i]));						
				Tf *= (1.0 - this->dx*this->temperature_slope); 
			}
		} else { 
			this->ODE.set_bc(i,Ti);
		}
	}
	printf("Tb=%lg E18 = %lg, energy = %lg\n",Tb,energy/(1e18*this->yb),energy*4.0*M_PI*pow(1e5*this->radius,2.0));
	return energy/(1e18*this->yb) - this->E18;
}

void Burst::calculate_cooling_curve(void) 
{
	int nsteps=0;
  	double er=0.0, en=0.0, eredd=0.0;   // eredd is the energy radiated with L<=L_edd

  	//double FEdd=3e10*this->g/(0.2*1.7);   // Eddington flux for solar
	double FEdd=3e10*this->g/0.2;   // Eddington flux for pure He
  	FEdd *= this->ZZ;   // we need to redshift this .. the idea is that any flux greater than the 
  	// Eddington flux *at infinity* is used to eject matter

  	double time_above_Edd=0.0;
  	double cooling_time=0.0;
	double last_time_output=0.0;
  	double dtstart=10.0;  // start timestep in sec

	printf("Running for time %lg seconds\n", this->time_to_run);
    
	// call the integrator
	clock_t timer;
	start_timing(&timer);
	printf("Starting integration\n");
	ODE.go(0.0, this->time_to_run, dtstart, 1e-8);
	stop_timing(&timer,"ODE.go");

	// extract lightcurve
	double *time = new double [ODE.kount];
	double *lum = new double [ODE.kount];
	for (int j=1; j<=ODE.kount; j++) {
		double flux=(this->g/env_g)*this->TEFF.get(ODE.get_y(1,j));
		time[j-1] = this->ZZ*ODE.get_x(j);
		lum[j-1] = 4.0*M_PI*pow(1e5*this->radius,2.0)*flux/(this->ZZ*this->ZZ);
	}
	this->lightcurve.minit(time,lum,ODE.kount);
	delete [] time;
	delete [] lum;
	
	// output	
	double lastFlux = 0.0;
	double maxeddrat = 0.0;
	if (this->output) {
		for (int j=1; j<=ODE.kount; j++) 
			output_result_for_step(j,FEdd,&lastFlux,&maxeddrat,&en,&er,&eredd,
				&cooling_time, &time_above_Edd,&last_time_output);
		fflush(this->fp); fflush(this->fp2);
	}
  	nsteps+=ODE.kount;
   	printf("number of steps = %d (total=%d) (time=%lg)  er=%lg eredd=%lg surface F/FEdd max=%lg\n", 
			ODE.kount, nsteps,this->time_to_run,er*4*M_PI*pow(this->radius*1e5,2.0), 
			eredd*4*M_PI*pow(this->radius*1e5,2.0), maxeddrat);
			
	if (this->output) {
	printf("Summary:\n");
	printf("E18=%lg, yb=%lg, yt=%lg, convect=%d\n", this->E18, this->yb, this->yt, this->convection_flag);
	printf("Energy input = %lg ergs\n", this->yb*this->E18*1e18*4.0*M_PI*pow(this->radius*1e5,2.0));
	printf("Total energy from surface = %lg,  as neutrinos=%lg\n", er*4.0*M_PI*pow(this->radius*1e5,2.0),en*4.0*M_PI*pow(this->radius*1e5,2.0));
	printf("Energy radiated at L<LEdd = %lg  (fraction of surface energy=%lg, effective E18=%lg)\n", eredd*4.0*M_PI*pow(this->radius*1e5,2.0), eredd/er, eredd/(this->yb*1e18));
	printf("where at infinity LEdd=%lg\n", FEdd*4.0*M_PI*pow(this->radius*1e5,2.0)/(this->ZZ*this->ZZ));

	printf("Observer sees (redshifted to infinity) energy %lg ergs\n", eredd*4.0*M_PI*pow(this->radius*1e5,2.0)/this->ZZ);
	printf("Flux dropped below FEdd after %lg s and time from FEdd to FEdd/e was %lg s\n", time_above_Edd,cooling_time);
	printf("Observer measures %lg and %lg s\n", time_above_Edd*this->ZZ, cooling_time*this->ZZ);
	}
	if (this->output)	
		fprintf(this->fp3, "%lg %lg %d %lg %d %lg %lg %lg %lg %lg\n", this->yb, this->E18, this->N, this->yt, this->convection_flag, 4.0*M_PI*pow(this->radius*1e5,2.0)*er, 
				4.0*M_PI*pow(this->radius*1e5,2.0)*en, 4.0*M_PI*pow(this->radius*1e5,2.0)*eredd, time_above_Edd, cooling_time);
}



void Burst::output_result_for_step(int j, double FEdd,
					double *lastFlux,double *maxeddrat,double *en,double *er,double *eredd,
					double *cooling_time, double *time_above_Edd,double *last_time_output)
{
	for (int i=1; i<=this->N+1; i++) calculate_vars(i,ODE.get_y(i,j),this->y[i],&this->CP[i],&this->K[i],&this->NU[i]);
	double T0;
	outer_boundary(ODE.get_y(1,j),this->K[1],this->CP[1],this->NU[1],&T0,&this->K[0],&this->CP[0],&this->NU[0]);
   	for (int i=2; i<=this->N+1; i++) this->F[i]=0.5*(this->K[i]+this->K[i-1])*(ODE.get_y(i,j)-ODE.get_y(i-1,j))/this->dx;
//this->F[1]=this->F[2];
	if (this->outer_boundary_flag) this->F[1]=(this->g/env_g)*this->TEFF.get(ODE.get_y(1,j));
	else this->F[1]=0.5*(this->K[1]+this->K[0])*(ODE.get_y(1,j)-T0)/this->dx;
   double lumn=0.0;
	for (int i=1; i<this->N; i++) lumn+=0.5*(this->y[i+1]-this->y[i])*(this->NU[i+1]+this->NU[i]);
	double dt;
	if (j==1) dt=ODE.get_x(j); else dt=ODE.get_x(j)-ODE.get_x(j-1);
	*er+=dt*this->F[1];
	if (this->F[1] < FEdd) *eredd+=dt*this->F[1]; else *eredd+=dt*FEdd;		
	*en+=lumn*dt;
	double FEdd_ph = 3e10*this->g*this->K[0]*this->y[0]/(3.03e-4*pow(T0,3.0));
	if (this->F[1]/FEdd_ph > *maxeddrat) *maxeddrat = this->F[1]/FEdd_ph;
	if (this->F[1]<FEdd && *lastFlux>FEdd) *time_above_Edd = ODE.get_x(j);
	*lastFlux = this->F[1];
	
	if (this->F[1]<FEdd/exp(1.0) && *time_above_Edd>0.0 && *cooling_time==0.0)
		*cooling_time = ODE.get_x(j) - *time_above_Edd;
	

	int n1=100;
	
	
	if (j==1 || (fabs(log10(fabs(this->ODE.get_x(j))*this->ZZ)-log10(fabs(*last_time_output))) >= 1e-2)) {

	// Note that in the following, time is now redshifted to infinity and the second column
	// is now Linfinity rather than the local flux at the star
   fprintf(this->fp2, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", 
		this->ZZ*ODE.get_x(j), 4.0*M_PI*pow(1e5*this->radius,2.0)*this->F[1]/(this->ZZ*this->ZZ), this->F[2], 
      	ODE.get_y(this->N-5,j), 0.0, this->F[10], *er, *en, lumn, *eredd, FEdd, FEdd_ph, ODE.get_y(1,j),
      	this->y[this->icool]);
      
	
	
 //  if (j % n1 == 0 || j==ODE.kount) {   // output every n1-th cycle
		// temperature profile
		fprintf(this->fp,"%lg\n",this->ZZ*ODE.get_x(j));
		double del;
		for (int i=1; i<=this->N; i++) {      
  			this->EOS.P=this->g*this->y[i]; this->EOS.T8=1e-8*ODE.get_y(i,j); this->EOS.rho=this->EOS.find_rho();
  			if (i>1) del = (ODE.get_y(i+1,j)-ODE.get_y(i-1,j))/(2.0*this->dx*ODE.get_y(i,j));
  			else del=1.0;
  			fprintf(this->fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", 
	  			this->y[i], ODE.get_y(i,j), this->F[i], this->NU[i],
	  			(this->F[i+1]-this->F[i])/(this->dx*this->y[i]), this->EOS.rho, this->EOS.CP()*this->EOS.rho, 
				ODE.get_d(i,j),1e8*pow(this->EOS.P/2.521967e17,0.25), this->EOS.opac(), del, this->EOS.del_ad(),
					2.521967e-15*pow(ODE.get_y(i,j),4)/(this->g*this->y[i]));	
	    }
		*last_time_output=(this->ODE.get_x(j))*this->ZZ;
	 }
}




double Burst::heat(double Tf)
{
	double energy = 0.0;
	double Ti = this->heat_Ti;
	this->EOS.P = this->g*this->heat_y;
	double dT = (Tf-Ti)*0.01;
	for (int i=1; i<=100; i++) {
		double T = Ti + dT*i;
		set_composition(this->shallow_composition_flag);
		this->EOS.T8 = T*1e-8;
		this->EOS.rho = this->EOS.find_rho();
		energy += dT * this->EOS.CP();	
	}
	return energy;	
}

void Burst::precalculate_vars(void) 
// calculate various quantities at each grid point as a function of temperature
// then during the run we can look them up in a table
{
	// the table is constructed in terms of log10(beta=Prad/P)
	this->nbeta=3000;
	this->betamin=-14.0;
	this->betamax=-0.001;
	this->deltabeta = (this->betamax-this->betamin)/(1.0*(this->nbeta-1));	

	this->rho_grid = matrix(this->N+2,this->nbeta);	
	this->CP_grid = matrix(this->N+2,this->nbeta);	
	this->K_grid = matrix(this->N+2,this->nbeta);	
	this->NU_grid = matrix(this->N+2,this->nbeta);	

	printf("Precalculating quantities...\n");

	for (int i=1; i<=this->N+1; i++) {
		
		this->EOS.P=this->g*this->y[i];
		if (i==icool) printf("Column for cooling = %lg\n",this->y[i]);	
		for (int j=1; j<=this->nbeta; j++) {		
			double beta = pow(10.0,this->betamin + (j-1)*(this->betamax-this->betamin)/(1.0*(this->nbeta-1)));
			if (this->y[i] > this->yb) set_composition(this->deep_composition_flag);
			else set_composition(this->shallow_composition_flag);
			this->EOS.T8 = 1e-8*pow(3.0*beta*this->EOS.P/7.5657e-15,0.25);
			this->EOS.rho = this->EOS.find_rho();
			this->rho_grid[i][j] = this->EOS.rho;
			this->CP_grid[i][j]=this->EOS.CP();
			double kappa=this->EOS.opac();
			//kappa = 1.0/((1.0/this->EOS.kappa_rad)+(1.0/this->EOS.potek_cond()));
		   this->K_grid[i][j]=3.03e20*pow(this->EOS.T8,3)/(kappa*this->y[i]);
			this->NU_grid[i][j]=this->EOS.eps_nu();	
			// add a cooling source in grid cell icool with luminosity L34*1e34 erg/s
			if (i==this->icool) this->NU_grid[i][j]+=1e34*this->L34*pow(this->EOS.T8/10.0,5.0)/(this->y[i]*this->dx*4.0*3.1415*pow(this->radius*1e6,2.0));
		}
	}
}


void Burst::set_composition(int composition_flag) {
	if (composition_flag==4) {
 		//printf("The cooling layer composition is 12C.\n");
 		this->EOS.A[1]=12.0; this->EOS.Z[1]=6.0;  // 28Si
	}	else
	if (composition_flag==3) {
		//printf("The cooling layer composition is 28Si.\n");
 		this->EOS.A[1]=28.0; this->EOS.Z[1]=14.0;  // 28Si
	} else
	if (composition_flag==2) {
		//printf("The cooling layer composition is Ca40.\n");
 		this->EOS.A[1]=40.0; this->EOS.Z[1]=20.0;  // 40Ca
	} else {
		//printf("The cooling layer composition is 56Ni.\n");
 		this->EOS.A[1]=56.0; this->EOS.Z[1]=28.0;  // 56Ni	
	}
  	this->EOS.set_Yi=this->EOS.Yi();  // hardwire the Yi's!  (increases speed)
  	this->EOS.set_Ye=this->EOS.Ye();
  	this->EOS.Z2=this->EOS.YZ2();
}


void Burst::set_up_grid(int ngrid)
{
	this->N=ngrid;   // number of grid points per decade in log column depth
	double ydeep=this->ydeep_factor*this->yb;// this is the target column depth of the innermost zone

  	// composition of the cooling layer
  	this->EOS.init(1);
	this->EOS.Q=900.0;
  	this->EOS.X[1]=1.0;   // pure something

  	// storage
  	this->CP=vector(this->N+2);
  	this->y=vector(this->N+2);
  	this->K=vector(this->N+2);
  	this->F=vector(this->N+2);
  	this->NU=vector(this->N+2);

 	// divide the entire region into this->N grid points
	int ib = (int) this->N * log(this->yb/this->yt)/log(ydeep/this->yt);
	this->dx=log(this->yb/this->yt)/(ib+0.5);    
	
  	for (int i=0; i<=this->N+1; i++) {
    	double x=log(this->yt)+this->dx*(i-1);
    	this->y[i]=exp(x);
  	}
  	printf("Grid has %d points (%d in the burning layer), delx=%lg\n", this->N, ib, this->dx);
}



void Burst::get_TbTeff_relation(double yHe, int comp)
// reads in the Flux-T relation from the data file 'grid_sorty'
// made by makegrid.cc
{
	double *temp, *flux;
	double y,T,F;
	// int npoints = 381;  // number of points to read in
	// temp = vector(npoints);
	// flux = vector(npoints);


	//  HOW TO MAKE THIS WORK?
	// if (this->envelope_file == 'none') {
	// 	char fname[50]="envelope_models/grid_sorty_";
	// 	if (comp == 0) {
	// 		sprintf(fname,"%sFe_%d",fname,(int) yHe);		
	// 	} else {
	// 		sprintf(fname,"%sSi_%d",fname,(int) yHe);
	// 	}
	// }	
	// else {
	// 	char fname = this->envelope_file;
	// }

	char fname[50]="envelope_models/grid_sorty_";
	if (comp == 0) {
		sprintf(fname,"%sFe_%d",fname,(int) yHe);		
	} else {
		sprintf(fname,"%sSi_%d",fname,(int) yHe);
	}
	printf("Reading TbTeff relation from file: %s\n",fname);


	FILE *fp = fopen(fname,"r");
	FILE *fp2;
	if (this->output) fp2=fopen("out/TbTeff", "w");

	// First run to find the number of points
	int npoints=0;
	while (!feof(fp)) {
		fscanf(fp, "%lg %lg %lg\n", &y, &T, &F);
		if (y == log10(this->yt)) ++npoints; // increment 1
	}
	fclose(fp);
	temp = vector(npoints);
	flux = vector(npoints);


	// Now reading the data
	fp = fopen(fname,"r");
	int count = 0;
	while (!feof(fp)) {
		
		fscanf(fp, "%lg %lg %lg\n", &y, &T, &F);
		if (y == log10(this->yt)) {  // select out the points which correspond to the top column
			count++; 
			temp[count] = pow(10.0,T);
			flux[count] = pow(10.0,F);
			// printf("%d %lg %lg %lg %lg %lg\n", count,y,T,F,temp[count],flux[count]);
			if (this->output) fprintf(fp2, "%d %lg %lg %lg %lg %lg\n", count,y,T,F,temp[count],flux[count]);
		}
	}
	printf("Read %d points from %s\n",count,fname);
		
	fclose(fp);
	if (this->output) fclose(fp2);
	
	this->TEFF.minit(temp,flux,count);
	
	free_vector(temp);
	free_vector(flux);
}



// --------------------------------------------------------------------------
// routines for the time integration

void Burst::derivs(double t, double T[], double dTdt[])
// calculates the time derivatives for the whole grid
{
	for (int j=1; j<=this->N; j++) calculate_vars(j,T[j], this->y[j], &this->CP[j], &this->K[j], &this->NU[j]);
  	outer_boundary(T[1],this->K[1],this->CP[1],this->NU[1],&T[0],&this->K[0],&this->CP[0],&this->NU[0]);
  	inner_boundary(T[this->N],this->K[this->N],this->CP[this->N],this->NU[this->N],&T[this->N+1],&this->K[this->N+1],&this->CP[this->N+1],&this->NU[this->N+1]);
  	for (int i=1; i<=this->N+1; i++)   this->F[i] = calculate_heat_flux(i,T);
  	for (int i=1; i<=this->N; i++) {
    	dTdt[i]=(this->F[i+1]-this->F[i])/(this->dx*this->CP[i]*this->y[i]);
		if (this->nuflag) dTdt[i]+=-(this->NU[i]/this->CP[i]);
  	}
  	dTdt[this->N+1]=0.0;
}


double Burst::calculate_heat_flux(int i, double *T)
{
	if (i==1 && this->outer_boundary_flag)
		return (this->g/env_g)*this->TEFF.get(T[i]);  //outer boundary flux
	else
		return 0.5*(this->K[i]+this->K[i-1])*(T[i]-T[i-1])/this->dx;	
}

double Burst::dTdt(int i, double *T)
// calculates the time derivative for grid cell i (used when calculating the jacobian)
{
  	int k=i-1; if (k<1) k=1;
  	int k2=i+1; if (k2>this->N+1) k2=this->N+1;
  	for (int j=k; j<=k2; j++) calculate_vars(j,T[j], this->y[j], &this->CP[j], &this->K[j], &this->NU[j]);
 	if (i==1) outer_boundary(T[1],this->K[1],this->CP[1],this->NU[1],&T[0],&this->K[0],&this->CP[0],&this->NU[0]);
 	if (i==this->N) inner_boundary(T[this->N],this->K[this->N],this->CP[this->N],this->NU[this->N],
		  								&T[this->N+1],&this->K[this->N+1],&this->CP[this->N+1],&this->NU[this->N+1]);
  	double f=(calculate_heat_flux(i+1,T)-calculate_heat_flux(i,T))/(this->dx*this->CP[i]*this->y[i]);
  	if (this->nuflag) f+=-(this->NU[i]/this->CP[i]);
  	return f;
}

void Burst::outer_boundary(double T1, double K1, double CP1, double NU1,
		double *T0, double *K0, double *CP0, double *NU0)  
{
	*T0=T1*(8.0-this->dx)/(8.0+this->dx);
	*K0=K1; *CP0=CP1;
	if (this->nuflag) *NU0=NU1;
}

void Burst::inner_boundary(double TN, double KN, double CPN, double NUN,
	double *TN1, double *KN1, double *CPN1, double *NUN1)
{
	*TN1=TN;
	*KN1=KN;
	*CPN1=CPN;  
	if (this->nuflag) *NUN1=NUN; else *NUN1=0.0;
}


void Burst::jacobn(double t, double *T, double *dfdt, double **dfdT, int n)
// calculates the Jacobian numerically
{
	double f, e=0.01;
  	// this assumes the arrays dfdt and dfdT are preinitialized to zero (I changed odeint to do this)
  	for (int i=1; i<n; i++) {
    	if (i>1) {
      	T[i-1]*=1.0+e; f=dTdt(i,T);
      	T[i-1]/=1.0+e; dfdT[i][i-1]=(f-dfdt[i])/(T[i-1]*e);
    	}
    	T[i]*=1.0+e; f=dTdt(i,T);
    	T[i]/=1.0+e; dfdT[i][i]=(f-dfdt[i])/(T[i]*e);
    	if (i<=n) {
      	T[i+1]*=1.0+e; f=dTdt(i,T);
      	T[i+1]/=1.0+e; dfdT[i][i+1]=(f-dfdt[i])/(T[i+1]*e);
    	}
  	}
}


void Burst::calculate_vars(int i, double T, double y, double *CP, double *K, double *NU)
{
	this->EOS.T8=1e-8*T; 
	this->EOS.P=this->g*y; 
	double beta = log10(2.521967e-15*pow(T,4)/this->EOS.P);
	
	// if beta lies outside the table, set it to the max or min value
	if (beta > this->betamax) beta = this->betamax;
	if (beta < this->betamin) beta = this->betamin;
		
	int j = 1 + (int) ((beta-this->betamin)/this->deltabeta);
	double interpfac=(beta-(this->betamin + (j-1)*this->deltabeta))/this->deltabeta;
	*K=this->K_grid[i][j] + (this->K_grid[i][j+1]-this->K_grid[i][j])*interpfac;
	this->EOS.rho=this->rho_grid[i][j] + (this->rho_grid[i][j+1]-this->rho_grid[i][j])*interpfac;
	*CP=this->CP_grid[i][j] + (this->CP_grid[i][j+1]-this->CP_grid[i][j])*interpfac;
	if (this->nuflag) *NU=this->NU_grid[i][j] + (this->NU_grid[i][j+1]-this->NU_grid[i][j])*interpfac; 
	else *NU=0.0;
}

