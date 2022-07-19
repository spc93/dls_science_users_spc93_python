import h5py
#compred to calc_coinc_1: remove -ve pixels (== bad pix mask); main loop as func; calc difference between two coincidence images
#roi should avoid dead pixels as -ve numbers look like bright spots - see module edges on image

#close('all')

datapath='/dls/i16/data/2014/cm4968-3/processing/'
#h5in=datapath+'twophoton_sum_b%i_%i.h5' % (450745, 450754)
h5in=datapath+'twophoton_sum_c458468_458601.h5'
infile = h5py.File(h5in,'r')
incorr2=infile['corr2']
insum=infile['sum']
ni, nj=195, 487

def circ(centre,radius,anglerad):
    #anglerad or radius can be vectors
    anglerad=array(anglerad)
    radius=array(radius)
    x=centre[0]+radius*sin(anglerad)
    y=centre[1]+radius*cos(anglerad)
    return array([x,y]).transpose().round()
def filled_square(centre, edge):
    return [[i,j] for i in range(centre[0]-edge/2, centre[0]+edge/2) for j in range(centre[1]-edge/2, centre[1]+edge/2)]



#ijvals=circ([70,350],20,pi*2*arange(0,1,.02)) #nothing
#ijvals=circ([70,350],50,pi*2*arange(0,1,.005)) #nothing
#mask1=filled_square([70,350],20); mask2=filled_square([130,350],20) #nothing
#mask1=filled_square([70,450],20); mask2=filled_square([130,450],20) #nothing
#mask1=filled_square([70,150],20); mask2=filled_square([130,150],20) #nothing
mask_list=[filled_square([ii,150],15) for ii in range(10,190,15)]

#add up images with these pixels from corr2 data
#load corr2 then
def conc_image(dat4d, mask):
    #dat4d is 4d coinc data set
    #mask is list if ij mask values
    ni, nj=195, 487
    imdat=zeros((ni,nj),'Int32')
    maskdat=zeros((ni,nj),'Int32')
    imsum=insum[...]
    nin=0;nout=0;sum_mask=0
    for ij in mask:
        print ij
        try:
            imdat+=incorr2[ij[0], ij[1]]
            #maskdat[ij[0], ij[1]]+=1
            #sum_mask+=incorr2[ij[0], ij[1], ij[0], ij[1]]
            nin+=1
        except:
            nout+=1
    print '%i pixels inside image, %i pixels outside image' % (nin, nout)
    imdat[imdat<0]=0 # remove -ve pixels
    return imdat

imdat0=conc_image(incorr2, mask_list[0])
imdiff_list=[]
for mask in mask_list[1:]:
    imdiff_list+=[conc_image(incorr2, mask)-imdat0]
    #figure(); pcolor(imdiff); axis('tight'); clim(-20, 20); colorbar() 
    
#plot and save
close('all')
ii=0; outfile='/home/spc93/tmp/coinc_3_frame%i.png'
for imdiff in imdiff_list:
    figure(); pcolor(imdiff); axis('tight'); clim(-20, 20); colorbar()
    savefig(outfile % ii)
    ii+=1
    
#imdat2=conc_image(incorr2, mask2)
#imdiff=imdat2-imdat1
#figure(); pcolor(imdat1); axis('tight'); colorbar() 
#figure(); pcolor(imdat2); axis('tight'); colorbar() 
#figure(); pcolor(imdiff); axis('tight'); clim(-20, 20); colorbar() 
#clim(-10, 10)