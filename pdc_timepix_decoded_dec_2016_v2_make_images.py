from dlstools.tictoc import tictoc
tt = tictoc()

#make clearer plots ##############################

#keep separate from main script
#run after main script



#%%
#get energy spectrum - total data for one file
#file0='long5_2016-12-12-15-56-24.091502_decoded.dat'
#file1='long5_2016-12-12-15-59-18.717693_decoded.dat'
#print tt
#file0 = allfiles['files0'][datadict[dataset][0]]    #    use just first file selected from file pairs
#file1 = allfiles['files1'][datadict[dataset][0]]

# takes over 200 sec so comment out
#t=pd.read_csv(datapath+datafolder_top+file0, sep='\t', usecols=[0,1,2,3], names=['I','J','ToA','ToT'])
#b=pd.read_csv(datapath+datafolder_bot+file1, sep='\t', usecols=[0,1,2,3], names=['I','J','ToA','ToT'])
#t.plot(y='ToT',kind='hist', bins=100)
#b.plot(y='ToT',kind='hist', bins=100)
#bad low-energy tail. top has worse energy resolution
#approx. chans for both: main 300-600; half 150-300
#top fwhm 360-540
#bot fwhm 340-440
#sum of two should be around 350-490

#%%
#analyse coincidence data from h5 
#### run first cell first


print tt
store = pd.HDFStore(datapath+'/'+dataset+'.h5')
all_coinc=store['all_coinc']; #load it back
isize=max(all_coinc.I+1)
jsize=max(all_coinc.J+1)


######### remove plotting for now....

#filter
roilist=sequential_roi_list(256,256,8,8)
print tt
for roi in roilist:
    #roi=roilist[36] #single roi in middle
    #time only
    #coinc=all_coinc[(all_coinc.I>=roi[0][0]) & (all_coinc.I<=roi[0][1]) & (all_coinc.J>=roi[1][0]) & (all_coinc.J<=roi[1][1]) ]
    #time and energy
    coinc=all_coinc[(all_coinc.I>=roi[0][0]) & (all_coinc.I<=roi[0][1]) & (all_coinc.J>=roi[1][0]) 
        & (all_coinc.J<=roi[1][1]) &  (all_coinc.ToT+all_coinc.ToT0>350)  & (all_coinc.ToT+all_coinc.ToT0<490)]
    #create image
    im_coinc=mk_image1(coinc, isize, jsize)
    plt.figure(); plt.imshow(im_coinc, vmax=2)

#coin.plot(y='ToT',kind='hist',bins=1+2*coinc_range, title='Pair No. '+str(filepairnum) ) #will crash if no coincidences
#coinc.plot(y='ToT',kind='hist', bins=100)
print tt

show()
