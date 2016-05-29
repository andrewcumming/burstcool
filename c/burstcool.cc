#include <stdio.h>
#include <math.h>
#include "../h/burst.h"

void calculate_chisq(Burst *burst,double distance);


int main(int argc, char *argv[])
{
	// Initialize the layer
	Burst burst;

	// turn output on or off
	burst.output = 1;	

	// Parse command line parameters
	burst.E18=atof(argv[1]);
	burst.yb=atof(argv[2]);
	burst.yt=atof(argv[3]);
	burst.convection_flag=atoi(argv[4]);
	burst.time_to_run=atof(argv[5]);
	if (argc > 6) burst.temperature_slope=atof(argv[6]);
	if (argc > 7) {
	//	burst.mass = atof(argv[7]);
	//	burst.radius = atof(argv[8]);
		burst.L34 = atof(argv[7]);
		burst.icool = atoi(argv[8]);
	}
	double distance=6.0;
	if (argc > 9) distance = atof(argv[9]);
	if (argc > 10) burst.output = atoi(argv[10]);
	
	// Setup the grid
	burst.setup();
	
	// Calculate the cooling
	burst.calculate_cooling_curve();

	// chi-squared
	calculate_chisq(&burst,distance);
}




void calculate_chisq(Burst *burst, double distance)
{
	// can add hard-coded data here
	// (time, luminosity and error in erg/s)
	//int nobs = 2;
	//double tobs[nobs] = { 100.0, 200.0 };
	//double Lobs[nobs] = { 1e38, 1e37 };
	//double eobs[nobs] = { 1e36, 1e35 };	
	#include "../data/1636_peak.cc"
//	#include "../data/1636.cc"

	double chisq=0.0;
	for (int i=0; i<nobs; i++) {
		chisq += pow((Lobs[i] - burst->lightcurve.get(tobs[i])/pow(distance/6.0,2.0))/eobs[i],2.0);
	}
	
	printf("chisq = %lf  (%lf)\n",chisq,chisq/nobs);	
}







