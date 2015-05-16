// A wrapper for the GSL root finder
//

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

double zbrent(double (*func)(double), double x1, double x2, double tol)
// Follows example at https://www.gnu.org/software/gsl/manual/html_node/Root-Finding-Examples.html#Root-Finding-Examples
{
	gsl_root_fsolver *s;
	int iter = 0, max_iter = 100, status;
	gsl_function F;
	// in the next line, cast the function to take an extra void* argument
	// as needed by GSL. 
	F.function = (double (*)(double,void*)) func;
	s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	gsl_root_fsolver_set (s, &F, x1, x2);
	
	do {
		iter++;
		status = gsl_root_fsolver_iterate (s);
		x1 = gsl_root_fsolver_x_lower (s);
		x2 = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x1, x2, 0.0, tol);
	} while (status == GSL_CONTINUE && iter < max_iter);
	
	gsl_root_fsolver_free (s);
	return gsl_root_fsolver_root(s);
}
