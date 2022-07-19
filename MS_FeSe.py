import sys
sys.path.append('/dls_sw/i16/software/python/DMS_Code')
sys.path.append('/dls_sw/apps/scisoftpy/2.7')
#===============================================================================
#         Gareth Nisbet Diamond Light Source - 25 June 2013
#===============================================================================
import os
import sys
import subprocess
from threading import Thread
import sys
from mpl_toolkits.mplot3d import Axes3D
dir = os.path.dirname(__file__)
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.dirname(__file__)+'calcms')
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib.pyplot import step, legend, xlim, ylim
from calcms import ts
import numpy as np
import scisoftpy.io as do
from scipy import ndimage
cif_file='FeSe_Cmma.cif'
hkllabels=0
hkl=[[3,3,0]]
hklint=hkl
azir=[1,0,0]
PV=[1,0]
energyrange=[7,7.2]
xmin, xmax =[-90,90]
numsteps=50
#===============================================================================
#                         DMS Calculation
#===============================================================================

mslist=[[np.NAN,np.NAN,np.NAN,np.NAN,np.NAN,np.NAN,np.NAN,np.NAN,np.NAN]] #full
################## Generate Reflist from Cif ###################################
SF, reflist, lattice , structure= ts.loadcif(cif_file,energyrange[-1])
refindex=~np.isnan(ts.vfind(reflist,np.round(hkl)-reflist).vindex())
SF=SF[refindex]
reflist=reflist[refindex]
SF=SF[ts.vfind(reflist,np.round(hkl)-reflist).vindex()]*SF
reflist=np.squeeze(reflist[np.where(SF!=0),:])
SF=np.squeeze(SF[np.where(SF!=0),:])
SF2=SF[ts.vfind(reflist,np.round(hkl)-reflist).vindex()]
loopnum=1

full, pv1, pv2, sfonly, pv1xsf1 = 0,0,0,1,0
#------------------------------------------------------------------------------ 
if pv1 + pv2 + sfonly + full >= 2:
    print('Choose only one intensity option')
    sys.exit()
elif pv1 + pv2 + sfonly + full + pv1xsf1 == 0:
    print('Geometry Only')
    mslist=[[np.NAN,np.NAN,np.NAN,np.NAN,np.NAN,np.NAN,np.NAN]]
for enval in np.linspace(energyrange[0],energyrange[1],numsteps):
    print(str(loopnum) + ' of ' +str(numsteps))

    #===========================================================================
    # SF0*Gauss*SF1*SF2*PV2
    #===========================================================================
    if full:
        calcstr='SF1*SF2*PV2'
        ms=ts.calcms(lattice,hkl,hklint,reflist,enval,azir,SF,SF2)#[:,[3,4,5]]
        polfull=ms.polfull(PV)
        mslist=np.concatenate((mslist,ms.polfull(PV)),0)
    #===========================================================================
    # PV1 only
    #===========================================================================
    elif pv1:
        calcstr='PV1'
        ms=ts.calcms(lattice,hkl,hklint,reflist,enval,azir,SF,SF2)
        mslist=np.concatenate((mslist,ms.pol1only(PV)),0)
    #===========================================================================
    # PV2 only
    #===========================================================================
    elif pv2:
        calcstr='PV2'
        ms=ts.calcms(lattice,hkl,hklint,reflist,enval,azir,SF,SF2)
        mslist=np.concatenate((mslist,ms.pol2only(PV)),0)
    #===========================================================================
    # SF only
    #===========================================================================
    elif sfonly:
        calcstr='SF1*SF2'
        ms=ts.calcms(lattice,hkl,hklint,reflist,enval,azir,SF,SF2)
        mslist=np.concatenate((mslist,ms.sfonly()),0)

        #===========================================================================
    # SF only
    #===========================================================================
    elif pv1xsf1:
        calcstr='SF1*PV1'
        ms=ts.calcms(lattice,hkl,hklint,reflist,enval,azir,SF)
        mslist=np.concatenate((mslist,ms.pv1xsf1(PV)),0)
    #===========================================================================
    # Geometry only - no structure factors
    #===========================================================================
    else:
        calcstr='Geometry Only'
        ms=ts.calcms(lattice,hkl,hklint,reflist,enval,azir)
        mslist=np.concatenate((mslist,ms.geometry()),0)

    loopnum=loopnum+1

keepindex=np.where([~np.isnan(mslist).any(1)])[1]
mslist=mslist[keepindex,:]

#===============================================================================
#                         MS Plotting
#===============================================================================
fig = plt.figure(figsize=(6, 3),dpi=130)
ax = fig.add_subplot(111,axisbg=[150/255.0,243/255.0,86/255.0])
mslist=np.array(mslist)
if pv1 + pv2 + sfonly + full + pv1xsf1 != 0:
    plt.scatter(mslist[:,3], mslist[:,7],c=mslist[:,-1],s=2,cmap=plt.cm.gray_r,lw = 0)
    plt.scatter(mslist[:,4], mslist[:,7],c=mslist[:,-1],s=2,cmap=plt.cm.gray_r,lw = 0)
    plt.colorbar()
else:
    plt.scatter(mslist[:,3], mslist[:,-1],s=2,lw = 0)
    plt.scatter(mslist[:,4], mslist[:,-1],s=2,lw = 0)
xlim(xmin,xmax)
ylim(energyrange[0],energyrange[1])
ax.set_xlim(xmin, xmax)
plt.title('CoCO3, hkl = ' + str(hkl)+' Incident polarization vector ' + str(str(PV)))
#plt.xlabel('Psi (deg)')
plt.xlabel(r'$\psi$ (deg)',fontsize=10)
plt.ylabel('Energy (keV)')
fig.subplots_adjust(bottom=0.2)
# plt.savefig('FeSe.svg',format='svg')
#plt.savefig('FeSe.pdf',format='pdf')
plt.show()
# ts.printfig('ps')
