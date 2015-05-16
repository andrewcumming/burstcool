#include <gsl/gsl_odeiv2.h>

class Ode_Int_Delegate {
public:
	virtual void derivs(double t, double T[], double dTdt[]){};
	virtual void jacobn(double, double *, double *, double **, int){};
};


class Ode_Int {
public:
	int ignore, kount, stiff, verbose, tri, use_gsl;
	double dxsav, minstep, hmax;
	void init(int n,Ode_Int_Delegate *delegate);
	void tidy(void);
	void go(double x1, double x2, double xstep, double eps);
	void go_gsl(double x1, double x2, long int nstep, double eps, int log_flag);
	void go_simple(double x1, double x2, int nstep);	
	void set_bc(int n, double num);
	double get_x(int i);
	double get_y(int n, int i);
	double get_d(int n, int i);
	double *xp, **yp;
	int nok, nbad;
	Ode_Int_Delegate *delegate;
	static int gsl_derivs (double t, const double y[], double dydt[], void * params);
	static int gsl_jacobn(double t, const double y[], double *dfdy,double dfdt[], void *params);
	
private:
	double *ynext,*derivs_y, *derivs_dydt, **derivs_dfdy;
	gsl_odeiv2_system sys;
	gsl_odeiv2_step *step;
	gsl_odeiv2_control *control;
	gsl_odeiv2_evolve *evolve;
	gsl_odeiv2_driver *driver; 

	double **dydxp,*hstr,*ystart;
	int kmax,nvar;
	void rkck(double y[], double dydx[], int n, double x, double h,
	   double yout[],
	   double yerr[]);
	void rkqs(double y[], double dydx[], int n, double *x, double htry, 
	   double eps,	double yscal[], double *hdid, double *hnext);
	void odeint(double ystart[], int nvar, double x1, double x2, double eps, 
	     double h1,double hmin, int *nok, int *nbad);
	void rk4(double y[], double dydx[], int n, double x, double h, double yout[]);
	void rkdumb(double vstart[], int nvar, double x1, double x2, int nstep);
	void rkscale(double vstart[], int nvar, double x1, double x2, double h1);

	double **d,*x;   // from stifbs.c

	void simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	    double xs, double htot, int nstep, double yout[]);
	void bansimpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	    double xs, double htot, int nstep, double yout[]);
	void trisimpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	    double xs, double htot, int nstep, double yout[]);
	void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	     double yscal[], double *hdid, double *hnext);
	void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);
	void lubksb(double **a, int n, int *indx, double b[]);
	void ludcmp(double **a, int n, int *indx, double *d);

	void tridag(double a[], double b[], double c[], double r[], double u[],
	     unsigned long n);
	void bandec(double **a, unsigned long n, int m1, int m2, double **al,
	     int *indx, double *d);
	void banbks(double **a, unsigned long n, int m1, int m2, double **al,
	     int *indx, double b[]);
};
