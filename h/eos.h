class Eos {
public:
	double potek_cond(void);
  // density, temperature and composition
  double rho, T8;
  // initialization
  void init(int n);
  void tidy(void);
  double *A, *Z, *X;
  // mean molecular weights
  double Ye(void);
  double Yi(void);
  double YZ2(void);
  double Z2, Yn;
  // equation of state
  double pe(void);
  double pemod(void);
  double ptot(void);
  double Chabrier_EF(void);
  double Fermi_Inv_1_2(double F);
  double FermiI(int k, double T8, double EF);
  // thermodynamics
  double f(void);
  double CP(void);
  double CV(void);
  double del_ad(void);
  double chi(double *x);
  double Gamma1(void);

  // opacity
  double eps_nu(void);
  double K_cond(double ef);
  double opac(void);
  double gff(double Z, double eta);
  double J(double x,double y);

  double lambda1, lambda2;
  double x(void);
  double gamma(void);
  double Uex(void);
  double eta(void);
  double lamei(int n);

  double TC(void);

  double Q, Q1,Q2,Q3,Q4, Qout, Q5;

  double cvion, cv_alpha, cvrad,cve;

  int accr;
  int gap;

  double Fep(int flag);
  double Eep(double q);

  double debug;
  double kes, kff, kcond, kgaunt, kappa_rad;

  double set_Ye, set_Yi;
  double RHO1,RHO2;

  double tmdrift(double dTdy);
  double Fermi(double n, double eta);
  void Fermi_derivs(double x, double ff[], double dfdx[]);
  static void Wrapper_Fermi_derivs(double x, double ff[], double dfdx[]);

  double find_rho(void);
  static double Wrapper_find_rho_eqn(double r);
  double find_rho_eqn(double r);
  double P;

  double s,w,L1,L2,I1,I2;
  double C12screen_factor, f_ee, f_ei;
  double f_eQ, f_ep;

  	double gamma_melt;

	double kncrit;

	double B;

 private:
  int ns;  // number of species
  double Fermi_n, Fermi_alpha;
	double expint(int n, double x);	

};
