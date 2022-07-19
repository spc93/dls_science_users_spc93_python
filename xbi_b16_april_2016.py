import dataloader
#matplotlib/pylab might need to be imported
t=dataloader.tiffloader('/media/sf_Dropbox/tmp/%05i.tif') #change the path to your own text
close('all')

crange=[100,500]# fix intensity range

for imdat in [2,3,5]:
    t(imdat); #lod image
    figure(); imshow(t.image0, cmap=gray()); clim(crange); colorbar(); axis('tight'); title(t.file)
    savefig('/media/sf_Dropbox/tmp/image%i.pdf' % imdat) #change the path to your own text

