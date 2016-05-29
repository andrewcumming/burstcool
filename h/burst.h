#include "../h/odeint.h"
#include "../h/spline.h"
#include "../h/eos.h"

class Burst: public Ode_Int_Delegate {
public:
	Burst();
	~Burst();

	void setup(void);
	void set_up_grid(int ngrid, int composition_flag);
	void calculate_cooling_curve(void);
	
	Spline lightcurve;
	
	int output;
  	double E18;  // energy input in units of 1e18 erg/g
	double yb, yt;
  	int convection_flag;
  	int nuflag;   // flag to switch on neutrino cooling
	double mass,radius;
	double temperature_slope;
	double time_to_run;

	void derivs(double t, double T[], double dTdt[]);
	void jacobn(double, double *, double *, double **, int);

	int icool;
	double L34;
	
private:
	FILE *fp,*fp2,*fp3;
	
	double dTdt(int i, double *T);
	void calculate_vars(int i, double T, double y, double *CP, double *K, double *NU);
	void outer_boundary(double T1, double K1, double CP1, double NU1, double *T0, double *K0, double *CP0, double *NU0);
	void inner_boundary(double TN, double KN, double CPN, double NUN, double *TN1, double *KN1, double *CPN1, double *NUN1);
	double calculate_heat_flux(int i, double *T);
	void precalculate_vars(void);
	void output_result_for_step(int j, double FEdd,
					double *lastFlux,double *maxeddrat,double*en,double*er,double*eredd,
					double *cooling_time, double *time_above_Edd,double *last_time_output);
	void get_TbTeff_relation(double yHe, int comp);

	double heat(double Tf);
	static double Wrapper_heat(double Tf);
	double set_initial_temperature_profile(double Tf);
	static double Wrapper_set_initial_temperature_profile(double Tf);
	
	Ode_Int ODE;
	Eos EOS;
	Spline TEFF;
	int N;   // grid cells
	double dx;  // grid
  	double *y, *CP, *K, *F, *NU;
	double **rho_grid, **CP_grid, **K_grid, **NU_grid;
	double betamin, betamax, deltabeta;
	double heat_Ti,heat_y;
	int nbeta;
  	double g, ZZ;
	int outer_boundary_flag;
  	int verbose;
};

