#include <stdio.h>
#include <math.h>

void set_ns_parameters(double mass, double radius, double *g, double *ZZ)
	// given mass in solar masses; radius in km, find gravity and redshift factor
{
	radius*=1e5;
	*ZZ=1.0/sqrt(1.0-2.0*6.67e-8*2e33*mass/(9e20*radius));
	*g=*ZZ*6.67e-8*2e33*mass/(radius*radius);	
	printf("NS parameters: M %lg M_sun, R %lg km, g14=%lg 1+z=%lg\n",
		mass, radius/1e5, *g/1e14, *ZZ);
}

void set_ns_redshift(double g, double R, double *mass, double *ZZ)
// given gravity g in cgs and radius in km
// sets the NS redshift and mass
{
	double y, x;
	R*=1e5;
	y = pow(R*g/9e20,2.0);
	// x is  GM/Rc^2
	x = y*(sqrt(1.0+y)-1.0);
	*mass = R*9e20*x/6.67e-8;
	*mass/=2e33;
	*ZZ=1.0/sqrt(1.0-2.0*x);	
}
