// class Spline
//
// Implements interpolation using the GSL cubic spline routines
// ( and also has a linear option)
//
// The data should be in the form of (x, y_i) rows, in binary 
// (double precision) with no seperators between the rows. 
// Spline::init reads x (assumed to be in column
// 1) and y (from column ycol). There are ncol columns and n rows.
//
// 'out_of_bounds_flag'. If the user sets this variable to zero,
// then 'get' returns zero if the x-coord is out of bounds; otherwise (default)
// it returns the first (or last) value in the table
//
// 'log_flag'. If this variable is set by the user to a non-zero
// value, the tabulated values are interpreted as log_10 values. Then 
// get returns 10^{table value}.
//
// getlin(x) does linear interpolation rather than spline
//
// size() returns the number of data points in the table
// (ie. it returns the value of this->num)
// get_x and get_y give access to the table
//

#include <stdio.h>
#include "math.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>

// definitions for fread & fwrite
#define DSIZE sizeof(double)
#define ISIZE sizeof(int)

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

double Spline::get_x(int i)
{
  return this->xtab[i];
}

double Spline::get_y(int i)
{
  return this->ytab[i];
}

int Spline::size(void)
{
  return this->num;
}

double Spline::get(double x)
{
  double u;

  if (x <= this->xtab[0]) {
    if (this->out_of_bounds_flag == 0) return 0.0;
    else {
      if (this->log_flag==0) return this->ytab[0];
      else return pow(10.0,this->ytab[0]); 
    }
  }

  if (x >= this->xtab[this->num-1]) {
    if (this->out_of_bounds_flag == 0) return 0.0;
    else {
      if (this->log_flag==0) return this->ytab[this->num-1];
      else return pow(10.0,this->ytab[this->num-1]); 
    }
  }

  u = gsl_spline_eval (spline, x, acc);

  if (this->log_flag==0) return u;
  else return pow(10.0, u);
}


void Spline::minit(double *x, double *y, int n)
  // initialize from memory rather than a file
{
  this->num=n;
  this->xtab=new double[n];
  this->ytab=new double[n];

  for(int i=0; i<n; i++) {  // copy x and y arrays
    this->xtab[i]=x[i+1];
    this->ytab[i]=y[i+1];
  }

  gsl_sort2(this->xtab,1,this->ytab,1,(size_t) n);
	
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init(spline,this->xtab,this->ytab,n);

  // default handling of out of bounds
  this->out_of_bounds_flag = 1;
  // default not logs
  this->log_flag=0;

  // 1st x-coordinate
  this->startx=this->xtab[0];
}

void Spline::tidy()
{
	delete [] this->xtab;
	delete [] this->ytab;
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
}





