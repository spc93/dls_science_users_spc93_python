from dlstools import dataloader
from dlstools.quickfit import *
from dlstools.dirty_fit import fit
import subprocess
path='/dls/i16/data/2015/cm12169-2/'
#create .dat file loader object
d=dataloader.dlsloader(path+'%i.dat')


#close('all') to close all plot windows
#load a .dat file...
#scans=range(515257,515283)
#for scan in scans:
#    print d(scan)
#    d.plot('eta','roi2_sum', hold=0); ylim(0,20000)
#    pause(1)
#    print d.sx
    
'''
d(515293)

subplot(2,2,1);d.plot('eta','roi2_sum', hold=1);
d.inc()
subplot(2,2,1);d.plot('eta','roi2_sum','r',hold=1);
d(515303)
subplot(2,2,2);d.plot('magrot','roi2_sum', hold=0);
d(515333)
subplot(2,2,3);d.plot('rotp','APD', hold=1);
d(515334)
subplot(2,2,3);d.plot('rotp','APD', hold=1);
savefig('/home/i16user/tmp/plot_dat.pdf')

#load the second image from it...
print p(2)
figure(); imshow(p.image_01); axis('tight'); title(p.file)

#load another .dat file...
d(477010)
#and now load the next one...
d.inc()
#quick way to plot with axis labels (see help: d.plot?)
d.plot('eta','sum')

#fit pseudo Voigt
fit(pv_c)
'''

import scisoftpy.io as do
def printfig():
    plt.savefig('/tmp/fig.ps',format='ps')
    cmdstr='lpr /tmp/fig.ps'
    subprocess.call(cmdstr,shell=True)
def orthogonal_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return numpy.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,-0.0001,zback]])
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
fig = plt.figure(1)
ax = fig.add_subplot(111,projection='3d')

magrot=0
v2=np.array([[]]*3).T

# for scan in range(515343,515343+72,2):
# scans=range(515343,515343+36,1)
# scans=range(515677,515713,1) # sigma-sigma
# scans=range(515714,515714+36,1) # sigma-pi
# scans=range(515751,515751+36,1)# total
# scans=range(515751,515785,1)# total
# scans=range(516017,516031,1)
# scans=range(516042,516042+36,1) # sigma-sigma
scans=range(516079,516079+36,1) # sigma - pi
scans=range(516116,516116+36,1) # total
# scans=range(515343+36,515343+36+36,1)
for scan in scans:
    d(scan)
#     v1=np.concatenate((np.array([np.ones(len(d.eta))*magrot]).T,np.array([d.eta]).T,np.array([d.APD]).T),1)
    v1=np.concatenate((np.array([np.ones(len(d.eta))*magrot]).T,np.array([d.eta]).T,np.array([d.roi1_sum]).T),1)
#     v1=np.concatenate((np.array([np.zeros(len(d.eta))+d.Ta[0]]).T,np.array([d.eta]).T,np.array([d.APD]).T),1)
#     v1=np.concatenate((np.array([d.xps3motor1]).T,np.array([d.eta]).T,np.array([d.APD]).T),1)
    v2=np.vstack([v2,v1])
    magrot=magrot+10
X = np.reshape(v2[:,0],(len(scans),len(d.eta)))
Y = np.reshape(v2[:,1],(len(scans),len(d.eta)))
Z = np.reshape(v2[:,2],(len(scans),len(d.eta)))
# dnp.plot.surface(Z,X[0,:],Y[0,:])
# ax.plot_surface(X, Y, Z,cmap=plt.cm.jet,rstride=1, cstride=1, linewidth=0, antialiased=True) 
ax.plot_surface(X, Y, Z,cmap=plt.cm.afmhot,rstride=1, cstride=1, linewidth=0, antialiased=True)
ax.view_init(45,45)
proj3d.persp_transformation = orthogonal_proj
# printfig()
plt.xlabel('eta (deg)')
plt.ylabel('APD')
# plt.title(d.metadata.cmd)
# plt.title('107 sigma/sigma @ Ta = 6K Energy = 5.223 ')
# plt.title('107 sigma/pi @ Ta = 6K Energy = 5.223 ')
plt.title('107 Total @ Ta = 6K Energy = 5.223 ')
'''
for scan in range(515993,515996,1):
    d=do.load('/dls/i16/data/2015/cm12169-2/'+str(scan)+'.dat')
    ax.plot(d.eta,d.APD)
axis('tight')
plt.xlabel('eta (deg)')
plt.ylabel('APD')
plt.title(d.metadata.cmd)
'''
