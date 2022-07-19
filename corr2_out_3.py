import h5py

datapath='/dls/i16/data/2014/cm4968-3/'
h5infmt=datapath+'processing/twophoton_%i.h5'
#h5infmt=datapath+'processing/twophoton_%i_%03i_%03i.h5'

goodscns=[458468, 458469, 458470, 458471, 458472, 458473, 458474, 458475, 458476, 458477, 458484, 458485, 458486, 458487, 458488, 458489, 458490, 458491, 458492, 458493, 458494, 458495, 458496, 458497, 458498, 458499, 458500, 458501, 458502, 458503, 458504, 458505, 458506, 458507, 458508, 458509, 458510, 458511, 458512, 458513, 458514, 458515, 458516, 458517, 458518, 458519, 458520, 458521, 458522, 458523, 458524, 458525, 458526, 458527, 458528, 458529, 458530, 458531, 458532, 458533, 458534, 458535, 458536, 458537, 458538, 458539, 458540, 458541, 458542, 458543, 458544, 458545, 458546, 458547, 458548, 458549, 458550, 458551, 458552, 458553, 458554, 458555, 458556, 458557, 458558, 458559, 458560, 458561, 458562, 458563, 458564, 458565, 458566, 458567, 458568, 458569, 458570, 458571, 458572, 458573, 458574, 458575, 458576, 458577, 458578, 458579, 458580, 458581, 458582, 458583, 458584, 458585, 458586, 458587, 458588, 458589, 458590, 458591, 458592, 458593, 458594, 458595, 458596, 458597, 458598, 458599, 458600, 458601]
scans=goodscns

ni, nj=195, 487

h5out=datapath+'processing/twophoton_sum_c%i_%i.h5' % (scans[0], scans[-1])
print '=== Summed data will be saved to '+ h5out
outfile = h5py.File(h5out,'w')
outcorr2 = outfile.create_dataset("corr2",(ni, nj, ni, nj), int16, chunks=(1, nj, ni, nj), compression='gzip', compression_opts=1) #two-event correlations
outsum = outfile.create_dataset("sum",(ni, nj), int32, compression='gzip', compression_opts=1) #two-event correlations

corr_2_chunk=0
sum_sum=0
goodscans=[]; badscans=[];
for scan in scans:
        try:
            h5in=h5infmt % scan
            infile = h5py.File(h5in,'r')
            incorr2=infile['corr2']
            insum=infile['sum']
            for ichunk in range(ni):
                print '=== Adding scan %i chunk %i to dataset' % (scan, ichunk)
                outcorr2[ichunk]+=incorr2[ichunk]
                outsum[...] += insum[...]
            infile.close()
            goodscans+=[scan]
        except:
            badscans+=[scan]

outfile.close()
print '=== Elvis has left the building'
print '=== Good scans: ', goodscans
print '=== Bad scans: ', badscans
print '=== Output file: %s' % h5out