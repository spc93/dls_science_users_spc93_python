import h5py
#December 2014
datapath='/dls/i16/data/2014/cm4968-5/'
h5infmt=datapath+'processing/twophoton_%i.h5'
#
goodscns=[490213, 490216, 490219, 490222, 490226, 490230, 490231, 490232, 490790, 490793, 490796, 490799, 490802, 490805, 490808, 490811, 490814, 490817, 490820, 490823, 490826, 490829, 490832, 490835, 490838, 490841, 490844, 490847, 490850, 490853, 490856, 490859, 490862, 490865, 490868, 490871, 490874, 490877, 490880, 490883, 490886, 490889, 490892, 490895, 490898, 490901, 490904, 490907, 490910, 490913, 490916, 490919, 490922, 490925, 490928, 490931, 490934, 490937, 490940, 490943, 490946]


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