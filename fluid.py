# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 21:31:31 2018

@author: Giacomo
"""

from fluid_functions import *

deb = 0;

## DATA
rotor = structtype();
rotor.p = 2; # number of pole pairs
mm = 1e-3; # millimeters
rotor.De = 200*mm; # [m], rotor outer diameter

rotor.Nb = 3; # number of flux-barriers
rotor.tb = np.array([4, 8, 15])*mm; # flux-barrier widths
rotor.wc = np.array([3, 7, 12, 10])*mm; # flux-carrier widths
rotor.Nstep = np.array([2, 4, 6]); # number of steps to draw the flux-barrier side
rotor.wrib_t = 1*mm; # [m], tangential iron rib width

# you can input flux-barrier angles or let the program compute them
#rotor.barrier_angles_el = np.array([14,26,38])*2; # [deg], electrical flux-barrier angles
#rotor.barrier_end = 'rect'; # choose 'rect' or comment

# you can define the rib width or comment
rotor.wrib = np.array([1,2,4])*mm; # [m], radial iron rib widths
# You can define the magnet width or comment
#rotor.wm = np.array([10,20,40])*mm;


## barrier points computation
barrier = structtype();
barrier = calc_fluid_barrier(rotor, deb);



## simple matlab plot
for bkk in range(0,np.size(barrier.X)):
    plt.plot(np.squeeze(barrier.X[bkk]), np.squeeze(barrier.Y[bkk]), '.-')
    
plt.axis('equal')
plt.show()


