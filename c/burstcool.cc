#include <stdio.h>
#include <math.h>
#include "../h/burst.h"

void calculate_chisq(Burst *burst);


int main(int argc, char *argv[])
{
	// Initialize the layer
	Burst burst;

	// Parse command line parameters
	burst.E18=atof(argv[1]);
	burst.yb=atof(argv[2]);
	burst.yt=atof(argv[3]);
	burst.ycomplete=atof(argv[4]);
	burst.convection_flag=atoi(argv[5]);
	burst.time_to_run=atof(argv[6]);
	burst.temperature_slope=atof(argv[7]);
	if (argc > 8) {
		burst.mass = atof(argv[8]);
		burst.radius = atof(argv[9]);
	}

	// turn output on or off
	burst.output = 1;	
	
	// Setup the grid
	burst.setup();
	
	// Calculate the cooling
	burst.calculate_cooling_curve();

	// chi-squared
	//calculate_chisq(&burst);
}




void calculate_chisq(Burst *burst)
{
	// can add hard-coded data here
	// 4U 1636 lightcurve
	int nobs = 1;
	double tobs[nobs] = { 100.0, 200.0 };
	double Lobs[nobs] = { 1e38, 1e37 };
	double eobs[nobs] = { 1e36, 1e35 };
		
	double chisq=0.0;
	for (int i=0; i<nobs; i++) {
		chisq += pow((Lobs[i] - burst->lightcurve.get(tobs[i]))/eobs[i],2.0);
	}
	
}







