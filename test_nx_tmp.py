import sys
sys.path.append('/dls_sw/i16/software/python')
from dlstools.pdnx import *
from matplotlib.pyplot import *



inpath = '/dls/i16/data/2020/mm23580-2/'
outpath = '/dls/science/users/spc93/misc_nexus_data/modified/'
filenum = 815893

p='/dls/science/users/spc93/misc_nexus_data/modified/%i.nxs'

n=pdnx(p % 815893)

#n.nx.entry1.scan.plot()

n.nx.plot()

root:NXroot
  @HDF5_Version = '1.10.4'
  @file_name = '/dls/science/users/spc93/misc_nexus_data/modif...'
  @file_time = '2020-04-01T12:57:04.454232'
  @h5py_version = '2.9.0'
  @nexusformat_version = '0.4.18'
  entry1:NXentry
    @default = 'scan'
    before_scan:NXcollection
      @target = '/entry1/before_scan'
      PPR:NXcollection
        ppchi = -44.999841517
        pppiezo1 = 0.0
        pppiezo2 = 12345.01236
        ppth1 = 0.0006406764
        ppth2 = -12.0003723606
        ppz1 = 11.0
        ppz2 = 11.0
      alpha:NXcollection
        alpha = -9.05796486487
      beamline_slits:NXcollection
        s1xcentre = -0.005
        s1xgap = 2.51
        s1ycentre = 0.0015
        s1ygap = 1.001
        s2xcentre = -0.00079375
        s2xgap = 4.9988625
        s2ycentre = 0.499
        s2ygap = 20.001
        s3xcentre = 2.69725
        s3xgap = 25.0015
        s3ycentre = -0.682
        s3ygap = 20.125
        s4xcentre = -5.875
        s4xgap = 36.6504
        s4ycentre = 0.0
        s4ygap = 28.005
        shtr3x = 11.93
        shtr3y = 4.25
      beta:NXcollection
        beta = 7.78667847136
      delta_offset:NXcollection
        delta_offset = 0.0
      dettrans:NXcollection
        dettrans = -27.0
      diffractometer_sample:NXcollection
        azih = 1.0
        azik = 0.0
        azil = 0.0
        beta = -10000000.0
        chi = 90.3482306912
        delta = 64.6931461735
        delta_axis_offset = 8.8
        en = 2.96699915623
        eta = 30.9881143025
        gam = 0.0
        h = -0.0282530656216
        k = -0.0189992729219
        kphi = 64.357
        l = 4.99932225101
        mu = 0.0
        phi = 6.76906790421
      dummypd:NXcollection
        x = 244.0
        y = 0.0
        z = 0.0
      gains_atten:NXcollection
        Atten = 0.0
        Transmission = 1.0
        diode_gain = 0
        ic1_gain = 2
        ic2_gain = 4
      jjslits:NXcollection
        s5xgap = 0.5
        s5xtrans = -0.0601
        s5ygap = 0.5002
        s5ytrans = 0.0
        s6xgap = 4.0786
        s6xtrans = 3.5176
        s6ygap = -0.5994
        s6ytrans = 8.5489
      lakeshore:NXcollection
        Ta = 47.099
        Tb = 47.625
        Tc = 3.15
        Td = 3.15
        Tset = 47.1
      mirrors:NXcollection
        m1piezo = 0.0
        m1pitch = 3.85035391186
        m1roll = 0.096548773231
        m1x = 1.6593395
        m1y = 5.2013125
        m1yaw = 0.399500886187
        m2bender = 305.0
        m2pitch = 3.85138940778
        m2roll = 0.399999978667
        m2x = 14.2
        m2y = 9.643875
        m2yaw = -0.037823791436
        m3pitch = 0.4497426
        m3x = -1.75e-05
        m4pitch = 0.4498818
        m4x = 2.1993175
      mono:NXcollection
        T1dcmSi111 = 116.840114372
        T2dcmSi111 = 11.7055480577
        bragg = -41.512175
        dcmfinepitch = 0.041053193124
        dcmlat = -41.4206
        dcmpitch = 25.5258322
        dcmroll1 = 96.2076116
        dcmroll2 = -22.9196032
        en = 2.96699915623
        perp = 40.79175
      mrwolf:NXcollection
        Day = 23
        Hours = 2
        Minutes = 19
        Month = 2
        Seconds = 21
        Year = 2020
      offsets:NXcollection
        base_z_offset = -13.8982188065
        idgap_offset = 0.529334663313
        m1y_offset = -15.7654858065
        m2_coating_offset = 11
        m2y_offset = -22.3229608065
        ztable_offset = -20.8213133065
      p2:NXcollection
        p2rot = 0.1845
        p2x1 = -5.02625
        p2y1 = 185.6515625
        p2zbot = 125.8384425
        p2ztop = -41.4271875
      pa:NXcollection
        stokes = 0.00012
        thp = -38.3771
        tthp = -76.199325
        zp = 4.8597
      positions:NXcollection
        Base_z = 7.06860475
        Base_z1 = 7.06860475
        Base_z2 = 7.06860775
        Base_z3 = 7.06860775
        base_y = 24.7925
        spara = 0.564166472585
        sperp = 1.36028507718
        sx = 0.3999
        sy = 1.4173
        sz = 8.97998
        table_horiz = 20.575
        table_vert = 0.14542
      ppchi:NXcollection
        ppchi = -44.999841517
      ppchitemp:NXcollection
        ppchitemp = 25.0
      pppitch:NXcollection
        pppitch = 0.0203
      ppth1:NXcollection
        ppth1 = 0.0006406764
      ppth1temp:NXcollection
        ppth1temp = 24.9
      ppth2:NXcollection
        ppth2 = -12.0003723606
      ppth2temp:NXcollection
        ppth2temp = 25.4
      ppx:NXcollection
        ppx = 10.2507
      ppy:NXcollection
        ppy = -0.7004
      ppyaw:NXcollection
        ppyaw = -0.0199
      ppz1:NXcollection
        ppz1 = 11.0
      ppz1temp:NXcollection
        ppz1temp = 24.7
      ppz2:NXcollection
        ppz2 = 11.0
      ppz2temp:NXcollection
        ppz2temp = 25.1
      psi:NXcollection
        psi = -80.0150097843
      s7xgap:NXcollection
        s7xgap = 0.5
      s7xtrans:NXcollection
        s7xtrans = -0.464
      s7ygap:NXcollection
        s7ygap = 0.4999
      s7ytrans:NXcollection
        s7ytrans = 0.0484
      source:NXcollection
        Uharmonic = 3
        idgap = 6.12395
        rc = 298.036590576
      xtlinfo:NXcollection
        UB11 = 1.11827489095
        UB12 = -0.327584532185
        UB13 = 0.00403231576852
        UB21 = 0.337165468331
        UB22 = 1.08686175454
        UB23 = -0.00177127478428
        UB31 = -0.0121600994705
        UB32 = 0.0100900058189
        UB33 = 0.321709932853
        a = 5.37915725585
        alpha1 = 90.0
        alpha2 = 90.0
        alpha3 = 90.0
        b = 5.53486462021
        c = 19.5287616747
    definition = 'NXmx'
    end_time = '2020-02-23T02:22:47.554Z'
      @target = '/entry1/end_time'
    entry_identifier = '815893'
    features = [                   3                    6 -4539647395532570642]
    instrument:NXinstrument
      BeamOK:NXpositioner
        beamOK = float64(81)
          @local_name = 'BeamOK.beamOK'
          @target = '/entry1/instrument/BeamOK/beamOK'
      atime:NXpositioner
        TimeSec = float64(81)
          @local_name = 'atime.TimeSec'
          @target = '/entry1/instrument/atime/TimeSec'
      atimetwo:NXpositioner
        TimeFromEpoch = float64(81)
          @local_name = 'atimetwo.TimeFromEpoch'
          @target = '/entry1/instrument/atimetwo/TimeFromEpoch'
      attenuator:NXattenuator
        attenuator_transmission = 1.0
      eta:NXpositioner
        eta = float64(81)
          @axis = '1'
          @label = '1'
          @local_name = 'eta.eta'
          @primary = '1'
          @target = '/entry1/instrument/eta/eta'
        soft_limit_max = 115.0
        soft_limit_min = -22.0
      ic1monitor:NXpositioner
        ic1monitor = float64(81)
          @local_name = 'ic1monitor.ic1monitor'
          @target = '/entry1/instrument/ic1monitor/ic1monitor'
      name = 'i16'
      pil3_100k:NXdetector
        calibration_date = '2018-07-11 11:15:52'
        calibration_scan_number = 705061
        count_time = float64(81)
          @local_name = 'pil3_100k.count_time'
          @target = '/entry1/instrument/pil3_100k/count_time'
        depends_on = '/entry1/instrument/pil3_100k/transformations/origin_offset'
        description = 'DetectorBase'
        id = 'detectorbase-1'
        image_data = [u'815893-pilatus3_100k-files/00001.tif', u'815893-pilat...
          @data_filename = 1
          @signal = 1
        maxval = float64(81)
          @local_name = 'pil3_100k.maxval'
          @target = '/entry1/instrument/pil3_100k/maxval'
        maxx = float64(81)
          @local_name = 'pil3_100k.maxx'
          @target = '/entry1/instrument/pil3_100k/maxx'
        maxy = float64(81)
          @local_name = 'pil3_100k.maxy'
          @target = '/entry1/instrument/pil3_100k/maxy'
        module:NXdetector_module
          data_origin = [0 0]
          data_size = [487 195]
          fast_pixel_direction = 0.172
            @depends_on = '/entry1/instrument/pil3_100k/module/module_offset'
            @offset = [0. 0. 0.]
            @transformation_type = 'translation'
            @units = 'mm'
            @vector = [ 0.70799834  0.00172762 -0.70621198]
          module_offset = 0.0
            @depends_on = '/entry1/instrument/pil3_100k/transformations/o...'
            @offset = [0. 0. 0.]
            @transformation_type = 'translation'
            @units = 'mm'
            @vector = [0. 0. 0.]
          slow_pixel_direction = 0.172
            @depends_on = '/entry1/instrument/pil3_100k/module/module_offset'
            @offset = [0. 0. 0.]
            @transformation_type = 'translation'
            @units = 'mm'
            @vector = [-2.1641950e-03  9.9999762e-01  2.7664700e-04]
        path = float64(81)
          @local_name = 'pil3_100k.path'
          @target = '/entry1/instrument/pil3_100k/path'
        saturation_value = 1000000
        sensor_material = 'Silicon'
        sensor_thickness = 0.32
        sum = float64(81)
          @local_name = 'pil3_100k.sum'
          @target = '/entry1/instrument/pil3_100k/sum'
        transformations:NXtransformations
          origin_offset = 1.0
            @depends_on = '/entry1/instrument/transformations/offsetdelta'
            @offset = [0. 0. 0.]
            @offset_units = 'mm'
            @transformation_type = 'translation'
            @units = 'mm'
            @vector = [ 62.21549545 -16.42424586 623.85767627]
        type = 'DetectorBase'
      rc:NXpositioner
        rc = float64(81)
          @axis = '1'
          @local_name = 'rc.rc'
          @target = '/entry1/instrument/rc/rc'
      roi2:NXdetector
        description = ''
        id = ''
        roi2_maxval = float64(81)
          @local_name = 'roi2.roi2_maxval'
          @target = '/entry1/instrument/roi2/roi2_maxval'
        roi2_maxx = float64(81)
          @local_name = 'roi2.roi2_maxx'
          @target = '/entry1/instrument/roi2/roi2_maxx'
        roi2_maxy = float64(81)
          @local_name = 'roi2.roi2_maxy'
          @target = '/entry1/instrument/roi2/roi2_maxy'
        roi2_sum = float64(81)
          @local_name = 'roi2.roi2_sum'
          @target = '/entry1/instrument/roi2/roi2_sum'
        type = ''
      source:NXsource
        current = 298.0366
          @units = 'mA'
        energy = 3.0009
          @units = 'GeV'
        name = 'DLS'
        probe = 'x-ray'
        type = 'Synchrotron X-Ray Source'
      transformations:NXtransformations
        delta = float64(81)
          @axis = '1'
          @depends_on = 'entry1/instrument/transformations/gamma'
          @local_name = 'kdelta.kdelta'
          @target = '/entry1/instrument/transformations/delta'
          @transformation_type = 'rotation'
          @units = 'deg'
          @vector = [0. 1. 0.]
        gamma = float64(81)
          @axis = '1'
          @depends_on = '.'
          @local_name = 'kgam.kgam'
          @target = '/entry1/instrument/transformations/gamma'
          @transformation_type = 'rotation'
          @units = 'deg'
          @vector = [1. 0. 0.]
        offsetdelta = float64(81)
          @axis = '1'
          @depends_on = 'entry1/instrument/transformations/delta'
          @local_name = 'delta_axis_offset.delta_axis_offset'
          @target = '/entry1/instrument/transformations/offsetdelta'
          @transformation_type = 'rotation'
          @units = 'deg'
          @vector = [ 0. -1.  0.]
    measurement:NXdata
      @axes = 'eta'
      @signal = 'roi2_sum'
      TimeFromEpoch = float64(81)
      TimeSec = float64(81)
      beamOK = float64(81)
      count_time = float64(81)
      delta_axis_offset = float64(81)
      eta = float64(81)
      ic1monitor = float64(81)
      kap = float64(81)
      kdelta = float64(81)
      kgam = float64(81)
      kmu = float64(81)
      kphi = float64(81)
      kth = float64(81)
      maxval = float64(81)
      maxx = float64(81)
      maxy = float64(81)
      path = float64(81)
      rc = float64(81)
      roi2_maxval = float64(81)
      roi2_maxx = float64(81)
      roi2_maxy = float64(81)
      roi2_sum = float64(81)
      sum = float64(81)
    pil3_100k:NXdata
      TimeFromEpoch -> /entry1/instrument/atimetwo/TimeFromEpoch
      TimeSec -> /entry1/instrument/atime/TimeSec
      beamOK -> /entry1/instrument/BeamOK/beamOK
      count_time -> /entry1/instrument/pil3_100k/count_time
      delta -> /entry1/instrument/transformations/delta
      eta -> /entry1/instrument/eta/eta
      gamma -> /entry1/instrument/transformations/gamma
      ic1monitor -> /entry1/instrument/ic1monitor/ic1monitor
      kappa -> /entry1/sample/transformations/kappa
      maxval -> /entry1/instrument/pil3_100k/maxval
      maxx -> /entry1/instrument/pil3_100k/maxx
      maxy -> /entry1/instrument/pil3_100k/maxy
      mu -> /entry1/sample/transformations/mu
      offsetdelta -> /entry1/instrument/transformations/offsetdelta
      path -> /entry1/instrument/pil3_100k/path
      phi -> /entry1/sample/transformations/phi
      rc -> /entry1/instrument/rc/rc
      sum -> /entry1/instrument/pil3_100k/sum
      theta -> /entry1/sample/transformations/theta
    program_name = 'GDA 9.15.0'
    roi2:NXdata
      TimeFromEpoch -> /entry1/instrument/atimetwo/TimeFromEpoch
      TimeSec -> /entry1/instrument/atime/TimeSec
      beamOK -> /entry1/instrument/BeamOK/beamOK
      delta -> /entry1/instrument/transformations/delta
      eta -> /entry1/instrument/eta/eta
      gamma -> /entry1/instrument/transformations/gamma
      ic1monitor -> /entry1/instrument/ic1monitor/ic1monitor
      kappa -> /entry1/sample/transformations/kappa
      mu -> /entry1/sample/transformations/mu
      offsetdelta -> /entry1/instrument/transformations/offsetdelta
      phi -> /entry1/sample/transformations/phi
      rc -> /entry1/instrument/rc/rc
      roi2_maxval -> /entry1/instrument/roi2/roi2_maxval
      roi2_maxx -> /entry1/instrument/roi2/roi2_maxx
      roi2_maxy -> /entry1/instrument/roi2/roi2_maxy
      roi2_sum -> /entry1/instrument/roi2/roi2_sum
      theta -> /entry1/sample/transformations/theta
    sample:NXsample
      beam:NXbeam
        incident_energy = 2.96700002543
          @axis = '1'
          @local_name = 'en.en'
        incident_wavelength = 0.417877288407
          @units = 'nm'
      depends_on = '/entry1/sample/transformations/phi'
      name = 'Default Sample'
      transformations:NXtransformations
        kappa = float64(81)
          @axis = '1'
          @depends_on = 'entry1/sample/transformations/theta'
          @local_name = 'kap.kap'
          @target = '/entry1/sample/transformations/kappa'
          @transformation_type = 'rotation'
          @units = 'deg'
          @vector = [ 0.          0.64278761 -0.76604443]
        mu = float64(81)
          @axis = '1'
          @depends_on = '.'
          @local_name = 'kmu.km'
          @target = '/entry1/sample/transformations/m'
          @transformation_type = 'rotation'
          @units = 'deg'
          @vector = [1. 0. 0.]
        phi = float64(81)
          @axis = '1'
          @depends_on = 'entry1/sample/transformations/kappa'
          @local_name = 'kphi.kphi'
          @target = '/entry1/sample/transformations/phi'
          @transformation_type = 'rotation'
          @units = 'deg'
          @vector = [0. 1. 0.]
        theta = float64(81)
          @axis = '1'
          @depends_on = 'entry1/sample/transformations/m'
          @local_name = 'kth.kth'
          @target = '/entry1/sample/transformations/theta'
          @transformation_type = 'rotation'
          @units = 'deg'
          @vector = [0. 1. 0.]
      ub_matrix = float64(1x3x3)
      unit_cell = float64(1x6)
        @angle_units = 'deg'
        @length_units = 'angstrom'
    scan:NXsubentry
      called_by = 'Not yet implemented in GDA'
      definition = 'NXclassicScan'
      end_time -> /entry1/end_time
      measurement:NXdata
        @axes = 'eta'
        @signal = 'roi2_sum'
        TimeFromEpoch = float64(81)
        TimeSec = float64(81)
        beamOK = float64(81)
        count_time = float64(81)
        delta_axis_offset = float64(81)
        eta = float64(81)
        ic1monitor = float64(81)
        kap = float64(81)
        kdelta = float64(81)
        kgam = float64(81)
        kmu = float64(81)
        kphi = float64(81)
        kth = float64(81)
        maxval = float64(81)
        maxx = float64(81)
        maxy = float64(81)
        path = float64(81)
        rc = float64(81)
        roi2_maxval = float64(81)
        roi2_maxx = float64(81)
        roi2_maxy = float64(81)
        roi2_sum = float64(81)
        scan_fields = [u'eta', u'beamOK', u'kphi', u'kap', u'kth', u'kmu', u'k...
        scan_header = [u' &SRS', u' SRSRUN=815893,SRSDAT=20200223,SRSTIM=02192...
        sum = float64(81)
      positioners:NXcollection -> /entry1/before_scan
        @target = '/entry1/before_scan'
        PPR:NXcollection
          ppchi = -44.999841517
          pppiezo1 = 0.0
          pppiezo2 = 12345.01236
          ppth1 = 0.0006406764
          ppth2 = -12.0003723606
          ppz1 = 11.0
          ppz2 = 11.0
        alpha:NXcollection
          alpha = -9.05796486487
        beamline_slits:NXcollection
          s1xcentre = -0.005
          s1xgap = 2.51
          s1ycentre = 0.0015
          s1ygap = 1.001
          s2xcentre = -0.00079375
          s2xgap = 4.9988625
          s2ycentre = 0.499
          s2ygap = 20.001
          s3xcentre = 2.69725
          s3xgap = 25.0015
          s3ycentre = -0.682
          s3ygap = 20.125
          s4xcentre = -5.875
          s4xgap = 36.6504
          s4ycentre = 0.0
          s4ygap = 28.005
          shtr3x = 11.93
          shtr3y = 4.25
        beta:NXcollection
          beta = 7.78667847136
        delta_offset:NXcollection
          delta_offset = 0.0
        dettrans:NXcollection
          dettrans = -27.0
        diffractometer_sample:NXcollection
          azih = 1.0
          azik = 0.0
          azil = 0.0
          beta = -10000000.0
          chi = 90.3482306912
          delta = 64.6931461735
          delta_axis_offset = 8.8
          en = 2.96699915623
          eta = 30.9881143025
          gam = 0.0
          h = -0.0282530656216
          k = -0.0189992729219
          kphi = 64.357
          l = 4.99932225101
          mu = 0.0
          phi = 6.76906790421
        dummypd:NXcollection
          x = 244.0
          y = 0.0
          z = 0.0
        gains_atten:NXcollection
          Atten = 0.0
          Transmission = 1.0
          diode_gain = 0
          ic1_gain = 2
          ic2_gain = 4
        jjslits:NXcollection
          s5xgap = 0.5
          s5xtrans = -0.0601
          s5ygap = 0.5002
          s5ytrans = 0.0
          s6xgap = 4.0786
          s6xtrans = 3.5176
          s6ygap = -0.5994
          s6ytrans = 8.5489
        lakeshore:NXcollection
          Ta = 47.099
          Tb = 47.625
          Tc = 3.15
          Td = 3.15
          Tset = 47.1
        mirrors:NXcollection
          m1piezo = 0.0
          m1pitch = 3.85035391186
          m1roll = 0.096548773231
          m1x = 1.6593395
          m1y = 5.2013125
          m1yaw = 0.399500886187
          m2bender = 305.0
          m2pitch = 3.85138940778
          m2roll = 0.399999978667
          m2x = 14.2
          m2y = 9.643875
          m2yaw = -0.037823791436
          m3pitch = 0.4497426
          m3x = -1.75e-05
          m4pitch = 0.4498818
          m4x = 2.1993175
        mono:NXcollection
          T1dcmSi111 = 116.840114372
          T2dcmSi111 = 11.7055480577
          bragg = -41.512175
          dcmfinepitch = 0.041053193124
          dcmlat = -41.4206
          dcmpitch = 25.5258322
          dcmroll1 = 96.2076116
          dcmroll2 = -22.9196032
          en = 2.96699915623
          perp = 40.79175
        mrwolf:NXcollection
          Day = 23
          Hours = 2
          Minutes = 19
          Month = 2
          Seconds = 21
          Year = 2020
        offsets:NXcollection
          base_z_offset = -13.8982188065
          idgap_offset = 0.529334663313
          m1y_offset = -15.7654858065
          m2_coating_offset = 11
          m2y_offset = -22.3229608065
          ztable_offset = -20.8213133065
        p2:NXcollection
          p2rot = 0.1845
          p2x1 = -5.02625
          p2y1 = 185.6515625
          p2zbot = 125.8384425
          p2ztop = -41.4271875
        pa:NXcollection
          stokes = 0.00012
          thp = -38.3771
          tthp = -76.199325
          zp = 4.8597
        positions:NXcollection
          Base_z = 7.06860475
          Base_z1 = 7.06860475
          Base_z2 = 7.06860775
          Base_z3 = 7.06860775
          base_y = 24.7925
          spara = 0.564166472585
          sperp = 1.36028507718
          sx = 0.3999
          sy = 1.4173
          sz = 8.97998
          table_horiz = 20.575
          table_vert = 0.14542
        ppchi:NXcollection
          ppchi = -44.999841517
        ppchitemp:NXcollection
          ppchitemp = 25.0
        pppitch:NXcollection
          pppitch = 0.0203
        ppth1:NXcollection
          ppth1 = 0.0006406764
        ppth1temp:NXcollection
          ppth1temp = 24.9
        ppth2:NXcollection
          ppth2 = -12.0003723606
        ppth2temp:NXcollection
          ppth2temp = 25.4
        ppx:NXcollection
          ppx = 10.2507
        ppy:NXcollection
          ppy = -0.7004
        ppyaw:NXcollection
          ppyaw = -0.0199
        ppz1:NXcollection
          ppz1 = 11.0
        ppz1temp:NXcollection
          ppz1temp = 24.7
        ppz2:NXcollection
          ppz2 = 11.0
        ppz2temp:NXcollection
          ppz2temp = 25.1
        psi:NXcollection
          psi = -80.0150097843
        s7xgap:NXcollection
          s7xgap = 0.5
        s7xtrans:NXcollection
          s7xtrans = -0.464
        s7ygap:NXcollection
          s7ygap = 0.4999
        s7ytrans:NXcollection
          s7ytrans = 0.0484
        source:NXcollection
          Uharmonic = 3
          idgap = 6.12395
          rc = 298.036590576
        xtlinfo:NXcollection
          UB11 = 1.11827489095
          UB12 = -0.327584532185
          UB13 = 0.00403231576852
          UB21 = 0.337165468331
          UB22 = 1.08686175454
          UB23 = -0.00177127478428
          UB31 = -0.0121600994705
          UB32 = 0.0100900058189
          UB33 = 0.321709932853
          a = 5.37915725585
          alpha1 = 90.0
          alpha2 = 90.0
          alpha3 = 90.0
          b = 5.53486462021
          c = 19.5287616747
      sample:NXsample
        name = 'Default Sample'
      scan_command -> /entry1/scan_command
      start_time -> /entry1/start_time
      title = 'test'
      user:NXuser -> /entry1/user01
        @target = '/entry1/user01'
        username = 'i16user'
      usergroup:NXcollection
        Bigears:NXuser
        Noddy:NXuser
        user01:NXuser
          username = 'i16user'
    scan_command = 'scan eta 30.9680320352948 32.5680320352948 0.02 BeamOK p...'
      @target = '/entry1/scan_command'
    scan_dimensions = 81
    scan_identifier = '26505c22-4dbe-44f7-be08-c6851e99cb05'
    start_time = '2020-02-23T02:19:17.647Z'
      @target = '/entry1/start_time'
    title = 'scan eta 30.9680320352948 32.5680320352948 0.02 BeamOK p...'
    title_old = 'Scan of sample with GDA'
    user01:NXuser
      @target = '/entry1/user01'
      username = 'i16user'
