#to do:
#

import sys
sys.path.append('/dls_sw/apps/scisoftpy/2.7')
sys.path.append('/dls_sw/i16/software/python')
from dlstools.dataloader import *
from dlstools import *
from matplotlib.pyplot import *

from dlstools.tictoc import tictoc
tt = tictoc()


out_folder = '/dls/science/users/spc93/pdc/'
close('all')

#%%
#Imports and set-up data used for all sections
#%matplotlib qt
#%matplotlib inline

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
pd.set_option('display.max_rows',18)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

min_entries=1000000 #files with fewer entries (photons) not used
datadict={'run1and5':range(2,106)+range(973,1009),
          '50M':range(4,108)+range(485,726),
          '100M':range(280,483)+range(720,972),
          'nobeam':range(112,276)} 
#datadict={'run1and5':range(22,25)} # small test set
#dataset='run1and5'
dataset='50M'

#datapath='/media/spc93/Data/data/I16_Timepix3/'
#datapath='/media/sf_data/I16_Timepix3/'
datapath='/dls/i16/data/2017/mt16227-1/processing/12_12_2016/Sorted'
#datafolder_top='/20161212_I16_W2J2_TOP_Decoded/'
#datafolder_bot='/20161212_I16_W5I6_BOT_Decoded/'
datafolder_top='/20161212_I16_W2J2_TOP_Decoded_Sorted/'
datafolder_bot='/20161212_I16_W5I6_BOT_Decoded_Sorted/'

def get_files_and_start_times(path,labelstr):
    
    #debug code -delete
    print path, labelstr
    ###################
    
    files=[]; times=[]
    file_list=os.listdir(path); file_list.sort();
    for file in file_list:
        firstline=pd.read_csv(path+file, sep='\t', usecols=[0,1,2,3], names=['I','J','ToA','ToT'], nrows=1)
        try:
            times+=[firstline.ToA[0]]
            files+=[file]
        except:
            print "no data: File:", file
    
    df=pd.DataFrame({'files'+labelstr:files, 'ToA'+labelstr:times})
    df['sec'+labelstr]=df['ToA'+labelstr]*2e-9;    df['mins'+labelstr]=df['sec'+labelstr]/60.0;    df['hours'+labelstr]=df['mins'+labelstr]/60.0;
    return df
   
def get_all_file_pairs_old(path_top, path_bot):      
    top_files_times=get_files_and_start_times(path_top,'0')
    bot_files_times=get_files_and_start_times(path_bot,'1')
    all=top_files_times.join(bot_files_times)
    all['ToA_offset']=all.ToA1-all.ToA0
    all['TimeOffsetSec']=all['ToA_offset']*2e-9
    return all


def get_all_file_pairs(path_top, path_bot):      
    top_files_times=get_files_and_start_times(path_top,'0')
    bot_files_times=get_files_and_start_times(path_bot,'1')
    bot_files_times_aligned=bot_files_times.copy()
    for ii in range(len(top_files_times)):
        sec0=top_files_times['sec0'][ii]
        idx=np.abs(bot_files_times.sec1-sec0).idxmin() #find index of closest start time of bot to top
        bot_files_times_aligned.loc[ii] = bot_files_times.loc[idx]
    all=top_files_times.join(bot_files_times_aligned)
    all['ToA_offset']=all.ToA1-all.ToA0
    all['TimeOffsetSec']=all['ToA_offset']*2e-9
    return all





#import numpy as np
#def find_nearest(array, value):
#    array = np.asarray(array)
#    idx = (np.abs(array - value)).argmin()
#    return array[idx]




def find_coincidences_old(path_top, path_bot, file_pair_table, file_num_list, coinc_time, coinc_range):
    #coinc_time=10   #count as coincidence if within this value
    #coinc_range=50  #plot range for time histograms
    num_events=0
    num_coinc=0
    for filepairnum in file_num_list:
        print "=== Starting file pair %i" % filepairnum
        t=pd.read_csv(path_top+file_pair_table.files0[filepairnum], sep='\t', usecols=[0,1,2,3], names=['I','J','ToA','ToT'])
        b=pd.read_csv(path_bot+file_pair_table.files1[filepairnum], sep='\t', usecols=[0,1,2,3], names=['I','J','ToA','ToT'])
        t['Module']=0; b['Module']=1
        b.I, b.J=255-b.I, 255*2-b.J
        tb=pd.concat([t,b]) #joint modules to make one big detector
        tb.sort_values('ToA',inplace=True)
        tb.ToA=tb.ToA-tb.ToA.min()
        tb['Tdiff']=tb.ToA-tb.ToA.shift(1)
        tb['Mdiff']=tb.Module-tb.Module.shift(1)
        tb['TdiffM']=tb.Tdiff/tb.Mdiff             #ToA difference: +ve if 0 first then 1; -ve if 1 then 0, inf if same module
        coinc_events_range=abs(tb.TdiffM)<=coinc_range #events within coinc_range (for plot)
        tb_coinc_range=tb[coinc_events_range]
        tb_coinc_range.plot(y='TdiffM',kind='hist',bins=1+2*coinc_range, title='Pair No. '+str(filepairnum) ) #will crash if no coincidences
        coinc_events=abs(tb.TdiffM)<=coinc_time #events within coinc_time (for coincidence data)
        tb_coinc=tb[coinc_events]
        tb_coinc_previous=tb.shift(1)[coinc_events]
        del tb_coinc_previous['Tdiff']; del tb_coinc_previous['Mdiff']; del tb_coinc_previous['TdiffM']; 
        tb_coinc_previous.rename(columns={'I':'I0','J':'J0','ToA':'ToA0','ToT':'ToT0','Module':'Module0'}, inplace=True)
        tb_both=tb_coinc.join(tb_coinc_previous)
        if num_events==0:
            all_coinc=tb_both.copy(deep=True)
        else:
            all_coinc=all_coinc.append(tb_both)
        num_events+=len(t)+len(b)
        num_coinc+=len(tb_both)
    all_coinc.reset_index()
    print "===Total events : %i\tCoincidence : %i\t Percentage of coincidences: %.3f" % (num_events, num_coinc, num_coinc*100./num_events)       
    return all_coinc.astype(int)



def find_coincidences(path_top, path_bot, file_pair_table, file_num_list, coinc_time, coinc_range, remove_adjacent=True, plot_histogram=False):
    #coinc_time=10   #count as coincidence if within this value
    #coinc_range=50  #plot range for time histograms
    #remove_adjacent: remove if there is a coincidence with the previous or next event (attempt at removing charge-sharing)
    #plot_histogram: plot coincidences as histogram (slow!)
    
    num_events=0
    num_coinc=0
    for filepairnum in file_num_list:
        print "=== Starting file pair %i" % filepairnum
        t=pd.read_csv(path_top+file_pair_table.files0[filepairnum], sep='\t', usecols=[0,1,2,3], names=['I','J','ToA','ToT'])
        b=pd.read_csv(path_bot+file_pair_table.files1[filepairnum], sep='\t', usecols=[0,1,2,3], names=['I','J','ToA','ToT'])
        t['Module']=0; b['Module']=1
        b.I, b.J=255-b.I, 255*2-b.J
        tb=pd.concat([t,b]) #joint modules to make one big detector
        tb.sort_values('ToA',inplace=True)
        tb.ToA=tb.ToA-tb.ToA.min()
        tb['Tdiff']=tb.ToA-tb.ToA.shift(1)
        tb['Mdiff']=tb.Module-tb.Module.shift(1)
        tb['TdiffM']=tb.Tdiff/tb.Mdiff             #ToA difference: +ve if 0 first then 1; -ve if 1 then 0, inf if same module
        tb['Prev_coinc']=abs(tb.ToA-tb.ToA.shift(2))<=coinc_range
        tb['Next_coinc']=abs(tb.ToA-tb.ToA.shift(-1))<=coinc_range
        coinc_events_range=abs(tb.TdiffM)<=coinc_range #events within coinc_range (for plot)
        tb_coinc_range=tb[coinc_events_range]
        if plot_histogram:
            tb_coinc_range.plot(y='TdiffM',kind='hist',bins=1+2*coinc_range, title='Pair No. '+str(filepairnum) ) #will crash if no coincidences
        if remove_adjacent: #coincidences are events between different modules where no coicidence in adjacent events
            #coinc_events=abs(tb.TdiffM)<=coinc_time & ~tb['Prev_coinc'] & ~tb['Next_coinc']
            coinc_events = (abs(tb.TdiffM)<=coinc_time) & (tb['Prev_coinc']==False) & (tb['Next_coinc']==False)
            tb['BeTrue'] = (abs(tb.TdiffM)<=coinc_time) & (tb['Prev_coinc']==False) & (tb['Next_coinc']==False)
        else: #coincidences are events between different modules 
            coinc_events=abs(tb.TdiffM)<=coinc_time #events within coinc_time (for coincidence data)
        
        tb_coinc=tb[coinc_events]
        print '=== No. coinc. events:', len(tb_coinc)
        tb_coinc_previous=tb.shift(1)[coinc_events]
        del tb_coinc_previous['Tdiff']; del tb_coinc_previous['Mdiff']; del tb_coinc_previous['TdiffM']; 
        del tb_coinc_previous['Prev_coinc']; del tb_coinc_previous['Next_coinc']; del tb_coinc_previous['BeTrue'];

        tb_coinc_previous.rename(columns={'I':'I0','J':'J0','ToA':'ToA0','ToT':'ToT0','Module':'Module0'}, inplace=True)

        print tb_coinc
        print tb_coinc_previous
        
        tb_both=tb_coinc.join(tb_coinc_previous)
        if num_events==0:
            all_coinc=tb_both.copy(deep=True)
        else:
            all_coinc=all_coinc.append(tb_both)
        num_events+=len(t)+len(b)
        num_coinc+=len(tb_both)
    all_coinc.reset_index()
    print "===Total events : %i\tCoincidence : %i\t Percentage of coincidences: %.3f" % (num_events, num_coinc, num_coinc*100./num_events)       
    return all_coinc.astype(int)
    #return tb ######## test



def sequential_roi_list(isize, jsize, idivisions, jdivisions):
    iwidth, jwidth=isize/idivisions, jsize/jdivisions 
    roilist=[]
    for ni in range(idivisions):
        for nj in range(jdivisions):  
            roilist+=[[[ni*iwidth,(ni+1)*iwidth-1],[nj*jwidth,(nj+1)*jwidth-1]]]
    return roilist
    
def mk_image1(event_table, isize, jsize):
    mtot=np.zeros([isize, jsize])
    for i in range(len(event_table)):
        mtot[event_table.I.iloc[i], event_table.J.iloc[i]]+=1   #event in roi
        mtot[event_table.I0.iloc[i], event_table.J0.iloc[i]]+=1 #event in coincidence
    return mtot

#%%
# get files corresponding to simultaneous pairs. Write to spreadsheet
print tt
allfiles=get_all_file_pairs(datapath+datafolder_top, datapath+datafolder_bot)
allfiles.plot('hours0','TimeOffsetSec',  marker='.' ); axis('tight')
plt.savefig(out_folder+'plot1.pdf')
writer=pd.ExcelWriter(out_folder+'allfiles.xlsx')
allfiles.to_excel(writer)
writer.save()
#%%
#create single h5 file for all coincidences for analysing
#all_coinc=find_coincidences(datapath+datafolder_top, datapath+datafolder_bot, allfiles, datadict[dataset], 10, 50)
print tt
all_coinc=find_coincidences(datapath+datafolder_top, datapath+datafolder_bot, allfiles, datadict[dataset][0:8], 10, 50, 
                            remove_adjacent=True, plot_histogram=False) #########test
store = pd.HDFStore(datapath+'/'+dataset+'.h5')
store['all_coinc']=all_coinc
store.close()
print tt
'''

#%%
#get energy spectrum - total data for one file
#file0='long5_2016-12-12-15-56-24.091502_decoded.dat'
#file1='long5_2016-12-12-15-59-18.717693_decoded.dat'
print tt
file0 = allfiles['files0'][datadict[dataset][0]]    #    use just first file selected from file pairs
file1 = allfiles['files1'][datadict[dataset][0]]

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
print tt
#coin.plot(y='ToT',kind='hist',bins=1+2*coinc_range, title='Pair No. '+str(filepairnum) ) #will crash if no coincidences
coinc.plot(y='ToT',kind='hist', bins=100)


show()


#%%


Notes:
energy spectrum of coincidences is broad, going up to ~500 whereas non-coinc. peaks at ~500


Future:
charge-sharing algorithm
remove multiple (>2) events as these are most likely cosmic rays

coco3 movie:
from matplotlib.animation import FuncAnimation
fig, ax = subplots(figsize=(5, 3))
axis('off');

cax = ax.imshow(images_e1[0], aspect=1/sin(29*pi/180), clim=(-50,800))
    
def animate(i):
    cax.set_array(images_e1[i])
    cax.set_array((images_e1[i]-images_e2[i])/images_e1[0])
    #ax.set_title('T=%.1f' % temp[i])
       
anim = FuncAnimation(fig, animate, frames=range(len(temp)))
#anim.save('day17sept17a_e1.mp4', dpi=400)

'''