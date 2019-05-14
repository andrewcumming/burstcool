Follows the thermal evolution of the outer layers of a neutron star after a thermonuclear flash. Based on the code used in Cumming & Macbeth (2004) and Cumming et al. (2006) ApJ papers.

To compile, first install

* GNU Scientific Library [GSL](http://www.gnu.org/software/gsl/) 
* [condegin13.f](http://www.ioffe.ru/astro/conduct/index.html) fortran routine by A. Potekhin to calculate thermal conductivity (put in directory `c`, you may need to remove deprecated 'pause' and 'stop' statements)

then

	mkdir o
	mkdir out
	make burstcool
	
should compile the code `burstcool` (you may need to change the compiler specified in the makefile to  whatever compiler you are using).

Parameters are given in the file `init.dat`:

* E18:	energy per gram in 10^18 erg/g   (roughly 1 MeV/nuc = 10^18 erg/g)
* yb:	base column depth (can be the actual number, or base 10 log, e.g. 1e12 or 12.0 will work)
* yt:	column depth at the top of the grid 
* burn:	a flag to indicate the initial temperature profile - 0=instantaneous burn (local deposition of energy), 1=adiabatic slope if the parameter <slope> is <0, or sets del=<slope> if <slope> is >0
* time_to_run:	seconds to run for (neutron star surface time)
* mass and radius: optional parameters to specify the mass and radius (you should give either both or none, default is 1.4 solar masses, 12 km)
* distance: units of kpc/10 km??
* output: Boolean for output to files (default 1)
* icool: Index of the grid cell for the cooling source (default 32)
* L34: Luminosity of the cooling source (default 0)
* ydeep_factor: default 100, y heating/base (heating at 1e12,base at 1e14 if ydeepfactor=100) (target column depth,
but based on grid construction its not exactly right)
* deep_composition : boolean
* shallow_composition : boolean

	
The code produces three output files in the directory `out`:

* `out/prof` -  one line per timestep giving luminosity etc., e.g. use this to plot the lightcurve
* `out/out`  -  full details of the layer structure as a function of time
* a line is added to `out/summary` with information such as the total energy radiated from the surface or in neutrinos etc.
* `info.txt` - full description of output files, variables and units
* More detailed information on all outputs is given in the `info.txt` file

`plot.pro` has IDL routines to make plots. To make a movie (uses ffmpeg):

	make cleanpng
	
	idl
	.com plot
	prof2, /png
	
	make movie      # not functionnal on current ffmpeg (?)
    
`plot.py` has python3 functions to make plots. To make a movie:

      make cleanpng
      python plot.py prof2
      make movie2

### Published lightcurves from this code

* [Cumming & Macbeth (2004)](http://lanl.arxiv.org/astro-ph/0401317) The Thermal Evolution following a Superburst on an Accreting Neutron Star
* [Cumming et al. (2006)](http://lanl.arxiv.org/astro-ph/0508432) Long Type I X-ray Bursts and Neutron Star Interior Physics
* [in 't Zand et al. (2014)](http://lanl.arxiv.org/abs/1312.5234) The Cooling Rate of Neutron Stars after Thermonuclear Shell Flashes
