#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../h/burst.h"

void calculate_chisq(Burst *burst,double distance);
void parse_parameters(char *fname,Burst &burst,double &);

int main(int argc, char *argv[])
{
	// Initialize the layer
	Burst burst;

	// turn output on or off
	burst.output = 1;	

	// Parse command line parameters
	// first determine the filename for the 'init.dat' parameter file
	char fname[200];
	char fnamedefault[10]="init.dat";
	switch(argc) {
		case 3:
			strcat(fname,"/tmp/init.dat.");
			strcat(fname,argv[1]);
			break;
		case 2:
			strcat(fname,"init/init.dat.");
			strcat(fname,argv[1]);
			break;
		default:
			strcat(fname,fnamedefault);
	}
	double distance=6.0;
	parse_parameters(fname,burst,distance);
	
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




void parse_parameters(char *fname,Burst &burst,double &distance) {
 	// Set parameters
	printf("============================================\n");
	printf("Reading input data from %s\n",fname);
	FILE *fp = fopen(fname,"r");
	char s1[100];
	char s[100];
	double x;				
	int commented=0;
	while (!feof(fp)) {   // we read the file line by line
		(void) fgets(s1,200,fp);		
		// ignoring lines that begin with \n (blank) or with # (comments)
		// or with $ (temperature profile)
		if (!strncmp(s1,"##",2)) commented = 1-commented;
		if (strncmp(s1,"#",1) && strncmp(s1,"\n",1) && strncmp(s1,">",1) && commented==0) {
			sscanf(s1,"%s\t%lg\n",s,&x);
			if (!strncmp(s,"E18",3)) burst.E18=x;
			if (!strncmp(s,"yb",2)) burst.yb=x;
			if (!strncmp(s,"yt",2)) burst.yt=x;
			if (!strncmp(s,"ydeep_factor",12)) burst.ydeep_factor=x;
			if (!strncmp(s,"burn",4)) burst.convection_flag=(int) x;
			if (!strncmp(s,"time_to_run",11)) burst.time_to_run=x;
			if (!strncmp(s,"slope",5)) burst.temperature_slope=x;
			if (!strncmp(s,"mass",4)) burst.mass=x;
			if (!strncmp(s,"radius",6)) burst.radius=x;
			if (!strncmp(s,"distance",8)) distance=x;
			if (!strncmp(s,"output",6)) burst.output=(int) x;
			if (!strncmp(s,"icool",5)) burst.icool=(int) x;
			if (!strncmp(s,"L34",3)) burst.L34=x;
			if (!strncmp(s,"deep_composition",16)) burst.deep_composition_flag=(int) x;
			if (!strncmp(s,"shallow_composition",19)) burst.shallow_composition_flag=(int) x;
			if (!strncmp(s,"L34",3)) burst.L34=x;
		}
	}

	fclose(fp);	
}
