// These functions to declare vectors and matrices
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

double *vector(long n)
{
	double *v = new double [n+1];
	return v;
}

double **matrix(long nr, long nc)
{
	double **m = new double*[nr+1];
	for (int i = 0; i < nr+1; ++i) m[i] = new double[nc+1];	
	return m;
}

void free_vector(double *v)
{
	delete [] v;
}

void free_matrix(double **m, long nr, long nc)
{
	for (int i = 0; i<nr+1; ++i) delete [] m[i];
	delete [] m;
}



