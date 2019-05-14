from numpy import array,where,pi
import matplotlib.pyplot as plt
import os 
import sys

def append_vars(line,varz,cols): # take line of file and append its values to variable lists 
    l=line.split()
    for var,col in zip(varz,cols):
        var.append(float(l[col]))

# Extract time and flux arrays
t,L,F = [],[],[]
with open('out/prof','r') as prof:
    for line in prof: append_vars(line,[t,L,F],[0,1,2])
t,L,F = array(t),array(L),array(F)

def lc():  # Plots the lightcurve 
    fig,ax = plt.subplots(1,1,figsize=(7.5,6))
    ax.set_xlabel(r'Time (s)',fontsize=14)
    ax.set_ylabel(r'Flux (erg cm$^{-2}$ s$^{-1}$)',fontsize=14)
    ax.loglog(t,F,'k-')
    if png: 
        plt.close()
        fig.savefig('lightcurve.png')
    else : plt.show()
        
def prof2(): # Makes movie frames of the temperature profile evolving over time '''

    if not os.path.exists('png'):   # Assuming code is being run from main directory 
        os.mkdir('png')
        
    # Plot window setup
    fig,[ax1,ax2] = plt.subplots(1,2,figsize=(15,6))
    ax1.set_xlabel(r'Time (s)',fontsize=14)
    ax1.set_ylabel(r'Flux (erg cm$^{-2}$ s$^{-1}$)',fontsize=14)
    ax2.set_xlabel(r'y (g cm$^{-2}$)',fontsize=14)
    ax2.set_ylabel(r'T (K)',fontsize=14) 
    if png: plt.close()
    ax1.loglog(t,F,'k-',lw=0.4,alpha=0.7)

    with open('out/out','r') as out:
        
        # Number of grid points
        ngrid = int(out.readline())
    #    print('Number of grid points =',ngrid)

        # Looping through datafile, plot frames
        count,time = 0,0
        y,T,FF,beta = [],[],[],[]
        for line in out:
            l = line.split()

            if len(l)==1:
                if count>0: 
                    
                    # plot current timestep
                    p = ax2.loglog(y,T,'k-')
                    fig.suptitle('t='+str(time))

                    # plot current location on lightcurve
                    t2,F2 = t[where(t<time)],F[where(t<time)]
                    ax1.loglog(t2,F2,'k-')
                    
                    if png:# save frame
                        fig.savefig(('png/%03d.png'%count))
                    else: plt.pause(0.0001)
                    
                    p.pop(0).remove()          # remove line
                    y,T,FF,beta = [],[],[],[]   # reset arrays at every timestep
                
                count += 1      
                time = float(l[0])
    #            print('t =',time)

            else: append_vars(line,[y,T,FF,beta],[0,1,2,12])

# Default is to save in png format instead of showing during exec.  Add '0' at end command line call to switch.
png = int(sys.argv[2]) if len(sys.argv)>2 else 1
if sys.argv[1] == 'lc' : lc()
elif sys.argv[1] == 'prof2': prof2()


def lc2(D=1):  # Plots the lightcurve with hour axis and observed fluence axis (give dist in kpc)
    fig,ax = plt.subplots(1,1,figsize=(10,6))
    ax.set_xlabel(r'Time (s)',fontsize=14)
    ax.set_ylabel(r'Luminosity (erg s$^{-1}$)',fontsize=14)
    ax.loglog(t,L,'k-')
    
    # Adding axes for hours and observed luminosity
    hours = [0.001,0.01,0.1,1,10,10]
    axtop = ax.twiny()
    axtop.loglog(t/3600,L,'b',alpha=0)
    axtop.set_xlabel(r'Time (hours)',fontsize=14) 
    axtop.set_xticks(hours)
    axtop.set_xticklabels(['{:g}'.format(hour) for hour in hours])

    kpc = 3.086e21 # cm
    axright = ax.twinx()
    axright.loglog(t,L/(4*pi*D**2*kpc**2),alpha=0)
    axright.set_ylabel(r'Observed fluence (erg cm$^{-2}$ s$^{-1}$)',fontsize=14) 

    if png: 
        plt.close()
        fig.savefig('lightcurve2.png')
    else : plt.show()

if sys.argv[1] == 'lc2' :
    if len(sys.argv)<4: sys.exit('Give distance to source in kpc')  
    lc2(D=float(sys.argv[3]))