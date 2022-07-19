
from dlstools import dataloader
path='/dls/i16/data/2018/cm19668-1/'
d=dataloader.dlsloader(path+'%i.dat')
p=dataloader.tiffloader(d, lambda obj: path+obj.pilatus2m_path_template)
close('all')
d(670844)

a=copy(d.sum)

sx=[]; sy=[]; r=[]
for sxval in arange(-2,2,.2):
    for syval in arange(-1,1,.2):
        sx+=[sxval]; sy+=[syval]; r+=[1./(sxval**2+syval**2+1)]


def vec2mat(vecx, vecy, vecz, n_inner=None):
    #matx, maty, matz = vec2mat(vecx, vecy, vecz, n_inner=None)
    #convert vectors from 2D scan to matrices
    #vecx,y,z: arrays (any dimension) or lists
    #matx,y,z: 2D arrays
    #n_inner: Number of points in inner loop - calculated if not specified
    #Arrays are truncated if the size doesn't match the required shape

    vx=array(vecx[:]); vy=array(vecy[:]); vz=array(vecz[:]) #get inputs in standard form
    if n_inner==None:   #calculate number in inner loop by looking for jumps
        jumps=abs(diff(vx)*diff(vy))
        n_inner=find(jumps>mean(jumps))[0]+1

    n_outer=len(vx)/n_inner
    #reshape matrices
    matx=vx[0:n_inner * n_outer].reshape(n_outer,n_inner)
    maty=vy[0:n_inner * n_outer].reshape(n_outer,n_inner)
    matz=vz[0:n_inner * n_outer].reshape(n_outer,n_inner)

    return matx, maty, matz

mx, my, mz = vec2mat(sx, sy, r)

figure()
#imshow(mz);
pcolor(mx, my, mz)


