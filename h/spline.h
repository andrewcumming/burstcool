#include <gsl/gsl_spline.h>

class Spline {
public:
  double get(double x);
  double get_x(int i);
  double get_y(int i);
  double startx;
  void minit(double *x, double *y, int n);
  void tidy();
  int out_of_bounds_flag;
  int log_flag;
  int size(void);
private:
  gsl_spline *spline;
  gsl_interp_accel *acc;
  double *ytab;
  double *derivs;
  double *xtab;
  int num;
};
