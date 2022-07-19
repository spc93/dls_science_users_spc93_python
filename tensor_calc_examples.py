#need this import if not in DAWN
from pylab import *
filestring=''
#spacegroup operator to have third (scalar) components +/- 1 for time reversal
#OK - now apply magnetic symmetry

# next: pass time symmetry keyword to symmetry functions
# symmeterize P=-1 T=-1 anapole in FeBO3 (hopefully it will exist and lie along a axis)
#functions to modify for time symmetry:
#apply_sym
#site_sym


#in calc_magSF_fourier(fourier_coeffs,sitevec,spacegroup,hkl) include extra parameter cells=[1,1,1] to represent the number of unit cells to use in each direction.
#need either a phase origin or to use complex coefficients
#use volumetric renderng programme to draw atoms, vectors and tensors in a unit cell.

#fix parity correction - done? check E1E2 K=1 term - gives errors
#determine if reflection is even, odd or mixed parity
#

from tools_for_tensor_calcs import *
set_printoptions(precision=3, suppress=True)
pol_theta=pi/4; #default tth_pol=90 deg unless specified
TimeEven=False #change to true for specific cases

plt.close('all')#close all existing plots



#CIFfile='/home/spc93/spc_cifs/GaFeO3_icsd_55840.cif'
#CIFfile='/home/spc93/Dropbox/spc_cifs/GaFeO3_icsd_55840.cif'
#CIFfile='/media/sf_Dropbox/spc_cifs/GaMnFeO3_174622.cif' # different setting - same as mcif? mag=a file editted to change gen pos to lowe case
# mag=c; elec= b; anapole=a
#hkl=array([0, 0, 3]); hkln=array([1,0,0]); lam=12.4/4.0; 
#Fe3=array([0.1525, 0.5827, 0.1893])
#lam=12.4/4.0; sitevec=Fe3; mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);# KK edge


#CIFfile='/home/spc93/Dropbox/spc_cifs/KOs2O6_icsd_162266.cif'
#hkl=array([4, 0, 2]); hkln=array([1,0,0]); lam=12.4/3.6;
#lam=12.4/5.0; #for Cs 
#sitevec=array([0.375, 0.375, 0.375]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);
#Os site
#sitevec=array([0., 0., 6.]); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/home/spc93/Dropbox/spc_cifs/BiFeO3_r3c_icsd_168319.cif'
#hkl=array([3, 3, 3]); hkln=array([1,0,0]); lam=12.4/7.1; 
#sitevec=array([0.2308, 0.2308, 0.2308]); mpol='E1E2';K=1; time=+1;parity=-1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/home/spc93/spc_cifs/KTP icsd_68705.cif'
#CIFfile='/media/sf_Dropbox/spc_cifs/KTP icsd_68705.cif'
#hkl=array([0, 0, 3]); hkln=array([1,0,0]); lam=12.4/4.0; 
#sitevec=array([.37807, .7806, .6880]); mpol='E1E1';K=2; time=+1;parity=1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);# KK edge
#lam=12.4/4.0; sitevec=array([0.37290, 0.5001, 0.99960]); mpol='E1E1';K=2; time=+1;parity=1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);# KK edge
#two K sites, two TYi sites


#CIFfile='/home/spc93/spc_cifs/ZrSiO4_icsd_100239.cif';
#hkl=array([0, 0, 14]); hkln=array([1,0,0]); lam=12.4/18.0; 
#sitevec=array([0, 3./4, 1./8]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/home/spc93/spc_cifs/KOs2O6_icsd_162266.cif'; # KOs2O6 Fd-3m pyrochlore
#hkl=array([0, 0, 2]); hkln=array([1,0,0]); lam=12.4/7.111 ; 
#sitevec=array([0, 0, 0]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG008=3.355/4; #pol_theta=arcsin(lam/2/d_an);
#sitevec=array([0, 0, 0]); mpol='E1E1';K=2; time=+1;parity=1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/home/spc93/spc_cifs/SmFeAsO_martinelli_icsd_163457.cif'; #SmFeAsO E1E2
#hkl=array([3, 2, 0]); hkln=array([0,0,1]); lam=12.4/7.111 ; 
#sitevec=array([3./4, 1./4, 1./2]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);
#sitevec=array([1./4, 1./4, .1372]);#Sm
#Sm1 Sm3+ 2 c 0.25 0.25 0.1372(1) 0. 1.
#O1 O2- 2 a 0.75 0.25 0. 0. 1.
#Fe1 Fe2+ 2 b 0.75 0.25 0.5 0. 1.
#As1 As3- 2 c 0.25 0.25 0.6599(2) 0. 1.

#CIFfile='/home/spc93/Dropbox/spc_cifs/Ho_icsd_56225.cif'; #Ho will not calculate cross section as not coded
#hkl=array([0, 0, 1]); hkln=array([1,0,0]); lam=12.4/9.659 ; sitevec=array([1./3, 2./3, 1./4]); mpol='E1E2';K=5; time=+1;parity=-1; ##d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/home/spc93/Dropbox/spc_cifs/ZnO Kisi et al icsd_67454.cif'; #ZnO E1E2
#CIFfile='/media/sf_Dropbox/spc_cifs/ZnO Kisi et al icsd_67454.cif'; TimeEven=True #add time symmetry (later)
#hkl=array([1, 1, 5]); hkln=array([0,0,1]); lam=12.4/9.659 ; sitevec=array([1./3, 2./3, 0.]); mpol='E1E2';K=3; time=+1;parity=-1; ##d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/media/sf_Dropbox/spc_cifs/HoFe2_icsd_103499.cif';   # HoFe2 E1E1
#CIFfile='/python/icsd_103499.cif';   # HoFe2 E1E1
#hkl=array([0, 0, 2]); hkln=array([1,0,0]);    lam=12.4/7.11; sitevec=array([0.625, 0.625, 0.625]); mpol='E1E1';K=2; time=+1;parity=+1; pol_theta=pi/4

#CIFfile='/python/Ge_icsd_53642.cif'; #Ge E1E2
#CIFfile='/media/sf_Dropbox/spc_cifs/Ge_icsd_53642.cif'; #Ge E1E2
#CIFfile='/home/spc93/Dropbox/spc_cifs/Ge_icsd_53642.cif'; #Ge E1E2
#TimeEven=True
#hkl=array([0, 0, 6]); hkln=array([1,0,0]); lam=12.4/11.103 ; sitevec=array([0., 0., 0.]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);
#hkl=array([0, 0, 2]); hkln=array([1,0,0]); lam=12.4/11.103 ; sitevec=array([0., 0., 0.]); mpol='E2E2';K=3; time=+1;parity=+1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/home/spc93/python/CuB2O4_icsd_20273.cif'; #CuB2O4 E1E1
#cu4b=array([0., 0., 0.5]); cu8d=array([0.0815, 0.25, 0.125]);
#hkl=array([0, 2, 2]); hkln=array([1,0,0]); lam=12.4/11.103 ; sitevec=cu8d; mpol='E1E1';K=2; time=+1;parity=1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/home/spc93/spc_cifs/Te_icsd_65692.cif'; #Te E1E1
#CIFfile='/home/spc93/python/Te_icsd_65692.cif'; #Te E1E1
#en_L1=4.94
#hkl=array([1, 0, 4]); hkln=array([1,0,0]); lam=12.4/en_L1 ; sitevec=array([-.2636, 0, 0.3333]); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/media/sf_Dropbox/spc_cifs/quartz_icsd_162490.cif'; #Quartz E1E1
#en=1.85
#hkl=array([0, 0, 1]); hkln=array([1,0,0]); lam=12.4/en ; sitevec=array([0.47000, 0, 2./3]); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/home/spc93/Dropbox/spc_cifs/berlinite_icsd_158614.cif'; #Berlinite E1E1
#en=1.57 #Al K
#hkl=array([0, 0, 1]); hkln=array([1,0,0]); lam=12.4/en ; sitevec=array([0.4160, 0.2590, 0.8840]); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/media/sf_Dropbox/spc_cifs/KDP_RT_I42d_icsd_173806.cif'; #KDP K(H2PO4) I-42d E1E2
#CIFfile='/media/sf_Dropbox/spc_cifs/KDP_ferro_Fdd2_icsd_68696.cif';#temp - parameters not checked
#en=2.145
#hkl=array([2, 2, 2]); hkln=array([1,0,0]); lam=12.4/en; sitevec=array([0.0, 0.0, 0.5]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG002=3.355; pol_theta=arcsin(lam/2/d_an);
#hkl=array([0, 0, 2]); hkln=array([1,0,0]); lam=12.4/en; sitevec=array([0.0, 0.0, 0.0]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG002=3.355; pol_theta=arcsin(lam/2/d_an);
#hkl=array([0, 0, 2]); hkln=array([1,0,0]); lam=12.4/en; sitevec=array([0.25, 0.25, 0.0125]); mpol='E1E1';K=2; time=+1;parity=1; d_an=d_PG002=3.355; pol_theta=arcsin(lam/2/d_an); # P site in low T ferro phase

#CIFfile='/home/spc93/spc_cifs/YIG_icsd_2012.cif'; #Fe E2E2
#CIFfile='/media/sf_Dropbox/spc_cifs/YIG_icsd_2012.cif';
#en_FeK=7.11; Fe16a=array([0, 0, 0]); Fe24d=array([3./8, 0, 1./4]); sitevec=Fe24d;
#hkl=array([8, 8, 8]); hkln=array([1,0,0]); lam=12.4/en_FeK ;  mpol='E2E2';K=4; time=+1;parity=+1; #d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/home/spc93/spc_cifs/GGG_sawada_icsd_84874.cif'; #Fe E2E2
#en_GdL1=8.38; Gd24c=array([1./8, 0, 1./4]); sitevec=Gd24c;
#hkl=array([0, 0, 8]); hkln=array([1,0,0]); lam=12.4/en_GdL1 ;  mpol='E2E2';K=4; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/home/spc93/Dropbox/spc_cifs/GaN_icsd_87830.cif';
#hkl=array([1, 1, 5]); hkln=array([0,0,1]); lam=12.4/10.367 ; sitevec=array([1./3, 2./3, 0.]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);
#hkl=array([-3, 0, 3]); hkln=array([1,0,0]); lam=12.4/10.367 ; sitevec=array([1./3, 2./3, 0.]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);

#use GaN CIF for CdSe
#CIFfile='/python/GaN_icsd_54698.cif';
#hkl=array([1, 1, 1]); hkln=array([0,0,1]); lam=12.4/3.538 ; sitevec=array([1./3, 2./3, 0.]); mpol='E1E2';K=3; time=+1;parity=-1; #d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);

#CdSe
#CIFfile='/home/spc93/spc_cifs/CdSe.cif';
#hkl=array([2, -1, 1]); hkln=array([0,0,1]); lam=12.4/3.538 ; sitevec=array([1./3, 2./3, 0.]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);

#FeBO3
#CIFfile='/media/sf_Dropbox/spc_cifs/FeBO3_icsd_34474.cif';
#CIFfile='/home/spc93/Dropbox/spc_cifs/FeBO3_icsd_34474.cif';
#hkl=array([-2,0,7]); hkln=array([1,0,0]); lam=12.4/5.2 ; sitevec=array([0., 0., 0.]); mpol='E2E2';K=4; time=+1;parity=+1; 
#lam=12.4/6.5 #Mn K
#d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);
#K=-1; magmode='explicit'; l_s=.5; sk=array([1,0,0]); lk=sk*l_s; e1k=sk;
#calculate magnetic structure factor by specifying fourier components of modulation
#K=-1; magmode='fourier'; tilt=.1*0; s=1; l=0; fourier_coeffs=[[array([0,0,-3]),cos(tilt)/2*array([1,0,0])],[array([0,0,0]),sin(tilt)*array([0,1,0])],[array([0,0,3]),cos(tilt)/2*array([1,0,0])]]; #tau in rlu and complex vector amplitude


#CoCO3
#CIFfile='/media/sf_Dropbox/spc_cifs/CoCO3_icsd_61066.cif'
#hkl=array([0, 0, 9]); hkln=array([1,0,0]); lam=12.4/5.2 ; sitevec=array([0., 0., 0.]); mpol='E2E2';K=4; time=+1;parity=+1; 
#lam=12.4/7.709 #Co K
#d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);
#K=-2 (E1); K=-1 nonres
#K=-1; magmode='explicit'; l_s=-2.0; sk=array([0,1,0]); lk=sk*l_s; e1k=sk;

#Nd2Ir2O7
#CIFfile='/media/sf_Dropbox/spc_cifs/Nd2Ir2O7_Fd3m_EntryWithCollCode156437.cif'
#CIFfile='/home/spc93/Dropbox/spc_cifs/Nd2Ir2O7_Fd3m_EntryWithCollCode156437.cif'
#sitevec=array([0., 0., 0.]); lam=12.4/11.215; #Ir site; Ir L3 edge
#sitevec=array([0.5, 0.5, 0.5]); lam=12.4/6.208; #Nd site; Nd L3 edge
#hkl=array([0, 0, 10]); hkln=array([1,0,0]);  sitevec=array([0., 0., 0.]); mpol='E1E1';K=2; time=+1;parity=+1; 
#TimeEven=True

#FeCrAs (right structure???)()
#CIFfile='/home/spc93/spc_cifs/FeCrAs_icsd_610175.cif'
#hkl=array([0.333, 0.333, 0.0]); hkln=array([1,0,0]); lam=12.4/5.2 ; sitevec=array([0.25, 0., 0.]); #Cr1
#K=-1; magmode='explicit'; l_s=0; sk=array([1,0,0]); lk=sk*l_s;
#d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/home/spc93/Dropbox/spc_cifs/tio2_rutile_icsd_64987.cif'; #Rutile
#en=3.61
#hkl=array([1, 2, 3]); hkln=array([1,0,0]); lam=12.4/en; sitevec=array([0.0, 0.0, 0.0]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG002=3.355; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/home/spc93/Dropbox/spc_cifs/CaIrO3_Cmcm_icsd_180015.cif'; #CaIrO3 Cmcm
#en=11.21
#hkl=array([0, 0, 5]); hkln=array([1,0,0]); lam=12.4/en; sitevec=array([0.0, 0.0, 0.0]); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG002=3.355; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/home/spc93/Dropbox/spc_cifs/Sr2IrO4_icsd_90741.cif'; #I41/acd
#en=11.21 #Ir L3
#hkl=array([1, 1, 28]); hkln=array([1,0,0]); lam=12.4/en ; sitevec=array([0.0, 1./4, 3./8]); mpol='E1E2';K=3; time=+1;parity=-1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/home/spc93/Dropbox/spc_cifs/La2CuO4_Pmab_icsd_155496.cif'; #I41/acd
#CIFfile='/home/spc93/Dropbox/spc_cifs/La2CuO4_Cmca_Yamada_icsd_65270.cif'; #Cmca
#en=8.979 #Cu K
#hkl=array([0, 0, 3]); hkln=array([1,1,0]); lam=12.4/en ; sitevec=array([0.0, 0.0, 0.0]); mpol='E1E1';K=2; time=+1;parity=1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/home/spc93/Dropbox/spc_cifs/Ca2Os2O7_icsd_97091.cif' #Imma #74
#en=10.87 #Os L3
#Os4c= [0.250, 0.250, 0.250]
#Os4b= [0, 0, 0.5]
#hkl=array([2, 1, 0]); hkln=array([1,1,0]); lam=12.4/en ; sitevec=array(Os4c); mpol='E1E1';K=2; time=+1;parity=1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/media/sf_Dropbox/spc_cifs/FeSe_Cmma.cif' #Low T orthorhombic
#TimeEven=True
#en=7.11 #FeK
#hkl=array([3, 3, 0]); hkln=array([1,1,0]); lam=12.4/en ; sitevec=array([0.25, 0, 0]); mpol='E2E2';K=4; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/media/sf_Dropbox/spc_cifs/CaBaCo4O7_icsd_162687.cif'#general positions changed to lower case in gedit
#TimeEven=True
#en=8 
#hkl=array([1,0,0]); hkln=array([1,1,0]); lam=12.4/en ; sitevec=array([0.0154, 0.9972, 0.9418]); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/media/sf_Dropbox/spc_cifs/Ca214_11K_Braden.cif'
#TimeEven=True
#en=2.838; #Ru L3 (L2=2.967)
#hkl=array([0,1,4]); hkln=array([0,1,0]); lam=12.4/en ; sitevec=array([0,0,0]); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/media/sf_Dropbox/spc_cifs/NaOsO3_scagnoli.cif'#general positions changed to lower case in gedit
#TimeEven=True
#en=10.8
#hkl=array([3,-2,0]); hkln=array([0,0,3]); 
#hkln=array([1.786,-0.015,0.285]);
#lam=12.4/en ; sitevec=array([0,0,0]); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/media/sf_Dropbox/spc_cifs/Cu2OSeO3.cif'
#TimeEven=True
#en=0.932 #CuL3
#site=Cu1= [0.8860, 0.8860, 0.8860] #4a
#site=Cu2=[0.1335, 0.1211, -.1281] #12b
#hkl=array([1,0,0]); hkln=array([0,0,1]); lam=12.4/en ; sitevec=array(site); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/media/sf_Dropbox/spc_cifs/ceosal.cif'
#TimeEven=True
#en=0.932 #CuL3
#Ce=[0.000000, 0.125700, 0.250000]
#Os=[0.250000, 0.250000, 0.000000] 
#hkl=array([0,1,0]); hkln=array([0,0,1]); lam=12.4/en ; sitevec=array(Os); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);


#FeAs
#CIFfile='/media/sf_Dropbox/spc_cifs/FeAs.cif';
#hkl=array([0, 0, 1]); hkln=array([1,0,0]); lam=12.4/5.2 ; 
#sitevec=array([0.0033, 0.1993, 0.25]); mpol='E1E2';K=3; time=+1;parity=-1; 
#lam=12.4/7.1#Fe K
#TimeEven=True
#d_an=d_PG008=3.355/4; pol_theta=arcsin(lam/2/d_an);
#K=-1; magmode='explicit'; l_s=.5; sk=array([1,0,0]); lk=sk*l_s; e1k=sk;
#calculate magnetic structure factor by specifying fourier components of modulation
#K=-1; magmode='fourier'; tilt=.1*0; s=1; l=0; fourier_coeffs=[[array([0,0,-3]),cos(tilt)/2*array([1,0,0])],[array([0,0,0]),sin(tilt)*array([0,1,0])],[array([0,0,3]),cos(tilt)/2*array([1,0,0])]]; #tau in rlu and complex vector amplitude

#CIFfile='/home/spc93/Dropbox/spc_cifs/GGG_sawada_icsd_84874.cif'
#TimeEven=True
#en=10.0 #ish
#Gd=[0.12500, 0.00000, 0.25000]
#hkl=array([0,2,0]); hkln=array([0,0,1]); lam=12.4/en ; sitevec=array(Gd); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/media/sf_Dropbox/spc_cifs/LuMnO3_lower_case_icsd_160377.cif'
#TimeEven=True
#en=1.5 #ish
#Lu=[/media/sf_Dropbox/spc_cifs0.51910, 0.58740, 0.25000]
#hkl=array([0,1,0]); hkln=array([0,0,1]); lam=12.4/en ; sitevec=array(Lu); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);

#LuMn)O3 orthorhombic phase
#CIFfile='/media/sf_Dropbox/spc_cifs/lumno3_orthorhombic.cif'
#TimeEven=True
#en=5.22 #ish
#Mn=[0.0, 0.0, 0.0]
#hkl=array([0,3.5,1]); hkln=array([0,0,1]); lam=12.4/en ; sitevec=array(Mn); 
#mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);
#K=-2 (E1); K=-1 nonres
#K=-1; magmode='explicit'; l_s=0.0; sk=array([0,1,0]); lk=sk*l_s; e1k=sk;

#CIFfile='/media/sf_Dropbox/spc_cifs/URu2Si2_cod_1538774_cifbib.cif'
#TimeEven=False
#en=5.22 #ish
#U=array([0,0,0]); Ru=array([0.5,0.25,0]);hkln=array([0,0,1]); lam=12.4/en ; 
#sitevec=U; 
#hkl=array([1,1,1])
#mpol='E1E1';K=2; time=-1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);
#K=-2 (E1); K=-1 nonres
#K=-1; magmode='explicit'; l_s=0.0; sk=array([0,1,0]); lk=sk*l_s; e1k=sk;

CIFfile='/home/spc93/spc_cifs/FeS2_Pa3_pyrites_cod_9000594_cifbib.cif'
TimeEven=True
en=2.47
hkl=array([0,0,1]); Fe=array([0,0,0]); S=array([0.385,0.385,0.385]);hkln=array([1,0,0]); lam=12.4/en ; 
sitevec=S; 
mpol='E1E1';K=2; time=+1;parity=1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);


#CIFfile='/media/sf_Dropbox/spc_cifs/Ca214_11K_Braden.cif'
#filestring='Ca214 013 Ru L3 edge'
#TimeEven=True
#en=2.838 #L3
##en=2.967 #L2
#hkl=array([0,1,3]); hkln=array([0,1,0]); lam=12.4/en ; sitevec=array([0,0,0]); mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);

#CIFfile='/media/sf_Dropbox/spc_cifs/nchem.2848-s6.cif'
#TimeEven=True
#en=8.0
#hkl=array([0,1,1]); Fe=array([0,0,0]); S=array([0.385,0.385,0.385]);hkln=array([0,2,0]); lam=12.4/en ; 
#sitevec=S; 
#mpol='E1E1';K=2; time=+1;parity=+1; d_an=d_PG004=3.355/2; pol_theta=arcsin(lam/2/d_an);



cif_obj=CifFile.CifFile(CIFfile); firstkey=cif_obj.keys()[0]; cb=cifblock=cif_obj[firstkey]; print cb.keys()
lattice=[float(cb['_cell_length_a'].partition('(')[0]),float(cb['_cell_length_b'].partition('(')[0]),float(cb['_cell_length_c'].partition('(')[0]),float(cb['_cell_angle_alpha'].partition('(')[0]), float(cb['_cell_angle_beta'].partition('(')[0]),float(cb['_cell_angle_gamma'].partition('(')[0])]
print 'Lattice:', lattice

try:
    symxyz=cb['_symmetry_equiv_pos_as_xyz']
except:
    symxyz=cb['_space_group_symop_operation_xyz'] #assume this is full group, not just generators
print symxyz


#spacegroup list (list of matrix/vector pairs in crystal  basis)
sglist=spacegroup_list_from_genpos_list(symxyz)

#force space group to be nonmagnetic by adding time reversed copy of each symmetry element
if TimeEven==True:
    print "=== Adding time-reversed symmetry operators'"
    sgnew=deepcopy(sglist)
    for sym in sgnew:
        sym[2]=-sym[2]
    sglist+=sgnew
###############################
#invert structure
#glist=invert_spacegroup(sglist)
#sitevec=-sitevec
################################


#calculate B matrix
B=latt2b(lattice)

n_sg=len(sglist); print 'Number of spacegroup operators: '+str(n_sg)
pglist= site_sym(sglist, sitevec); n_pg=len(pglist)
print 'Number of symmetry operators at atomic site: '+str(n_pg)
print 'Number of equivalent sites in unit cell: '+str(n_sg*1.0/n_pg)
print 'Equivalent atomic sites:'
print equiv_sites(sglist, sitevec)
crystalpglist= crystal_point_sym(sglist)
print 'Number of symmetry operators for point group of crystal: '+str(len(crystalpglist))

theta=arcsin(lam*norm(dot(B, hkl))/2)  #for diffraction
#theta=0.0;    #for transmission with hkl perp to E and q
print 'theta(deg)=', theta*180/pi

sym_phases=SF_symmetry(sitevec, hkl, sglist)            #return parity odd/even etc and display messages

##### experimental - comment or remove this section ##### using now for hofe2
#magvec=array([0,0,1])
#Tsym=T_symmetry_consistent_with_fm_vector(sglist,magvec)
#Tsym=T_symmetry_favouring_fm_vector(sglist,magvec)
#Tsym=ones(len(Tsym)) #all +ve
#sglist=intersection(sglist,gen_twelve_fold_about_v([1,0,0])) #new spacegroup - intesection of sg with field
#mag_sg=set_spacegroup_timesym(sglist,Tsym)
#re-do spacegroup and point group with new time signatures
#sglist=set_spacegroup_timesym(sglist,Tsym)
#sym_phases=SF_symmetry(sitevec, hkl, sglist)# info on new low-symmetry group

if not isGroup(sglist): print '<><><> This is NOT a group!!!'
#crystalpglist= crystal_point_sym(sglist)
#pglist= site_sym(sglist, sitevec);
#print 'Magnetic signatures: ', Tsym

##############################################




if K>=0:
    print 'Tensor calculation'
    Ts=list(rand(2*K+1));     #use random numbers for quick and dirty solution... (replace each element with 2K+1 element list in future)
    #Ts=[0.123456+0.392957*1.J, 0.836836+0.296849*1.J, 1, -0.836836+0.296849*1.J,  0.123456-0.392957*1.J]; #fixed pseudo random K=2 (Te etc)

    #############
    for Qp in range(1,K+1): #next four optional lines modify previous line to adopt Brouder's relationship between tensor components
        Qn=-Qp; rndr=rand(); rndi=rand();
        Ts[K-Qp]=rndr+rndi*1.J;
        Ts[K+Qp]=(-1)**Qp*(rndr-rndi*1.J);
    #############

    Tc1=spherical_to_cart_tensor(Ts)   #convert to cartesian tensor of same rank
    #Tc1=rand(3,3); print 'XXXXXX remove this line! generate random cartesian K=2 tensor'
    Tc_atom=apply_sym(Tc1, pglist, B, P=parity, T=time);  #apply site symmetry using site point group ops and B matrix
    Tc_crystal=apply_sym(Tc1, crystalpglist, B, P=parity, T=time);  #apply site symmetry using crystal point group ops and B matrix
    Fc=norm_array(calc_SF(Tc1, sitevec, hkl, sglist, B, P=parity, T=time));   #calc SF Crt tensor using crystal space group and B matrix
    
    #######for diagnostics
    #for sym_phase in sym_phases[0]:
    #   mat=sym_phase[0]; phases=sym_phase[1];
    #    #print mat; print  mean(phases); print dot(mat, mat); print;
    #   print mat; print phases; print;
    #print sym_phases[0]
    ##################
    Ts_atom=norm_array(cart_to_spherical_tensor(Tc_atom));    #atomic spherical tensor
    print 'Atomic spherical tensor:'; print Ts_atom
    Ts_crystal=norm_array(cart_to_spherical_tensor(Tc_crystal));    #crystal spherical tensor
    print 'Bulk crystal spherical tensor:'; print Ts_crystal
    Fs=norm_array(cart_to_spherical_tensor(Fc));    #SF spherical tensor
    print 'Structure factor spherical tensor:'; print Fs
elif magmode=='fourier':
    #calculate magnetic SF from fourier components
    mk=calc_magSF_fourier(fourier_coeffs,sitevec,sglist,hkl)
    if K==-1:
        #non-res mag using fourier
        sk=s*mk; lk=l*mk;
        print "=== Non-resonant magnetic diffraction"; print "|Fs|^2 = ",norm(sk)**2
    elif K==-2:
        #res mag using fourier
        e1k=mk;
        print "=== Resonant (E1E1) magnetic diffraction"; print "|Mk|^2 = ",norm(e1k)**2

#unit vetors in diffraction (theta) coordinate system
h=x=array([1,0,0])
q0=array([-sin(theta), cos(theta),0])
q1=array([sin(theta),cos(theta),0])
esig=z=array([0,0,1])
epi0=array([cos(theta),sin(theta),0])
epi1=array([cos(theta),-sin(theta),0])
y=array([0,1,0])

#Stokes parameters (incient beam) P3=lin, P1=45, P2=circ
P3=1; P1=P2=0;
mu=1./2*array([[1.+P3, P1-1.J*P2], [ P1+1.J*P2, 1.-P3]])

psilist=[]; totlist=[]; p0list=[]; p45list=[]; p90list=[]; p135list=[]; circplist=[]; circmlist=[]; Isslist=[]; Isplist=[]; Iabslist=[]; Itotabs=[]; f0_r=[]; f0_im=[]; Ipslist=[]; Ipplist=[]; Islist=[]; Iplist=[];
for psideg in range(361):
    psi=psideg*pi/180
    Uctheta=theta_to_cartesian(hkl,hkln,psi,B)
    #unit vetors in crystal Cartesian coordinate system
    h_c=x_c=dot(Uctheta, h)
    q0_c=dot(Uctheta,q0)
    q1_c=dot(Uctheta,q1)
    esig_c=z_c=dot(Uctheta,esig)
    epi0_c=dot(Uctheta,epi0)
    epi1_c=dot(Uctheta,epi1)
    y_c=dot(Uctheta,y)
    e45_c=(esig_c+epi1_c)/sqrt(2)
    e135_c=(esig_c-epi1_c)/sqrt(2)
 
    if K>=0:
        #scattering matrix G, analyser matrix A
        X_ss=Xtensor(mpol, K, time, parity, esig_c, esig_c,  q0_c,  q1_c);    f_ss=scalar_contract(X_ss, Fs);
        X_sp=Xtensor(mpol, K, time, parity, esig_c, epi1_c,  q0_c,  q1_c);    f_sp=scalar_contract(X_sp, Fs);
        X_ps=Xtensor(mpol, K, time, parity, epi0_c, esig_c,  q0_c,  q1_c);    f_ps=scalar_contract(X_ps, Fs);    
        X_pp=Xtensor(mpol, K, time, parity, epi0_c, epi1_c,  q0_c,  q1_c);    f_pp=scalar_contract(X_pp, Fs);
        #transmitted beam following incident beam and scattered beam

        X_ss_q0=Xtensor(mpol, K, time, parity, esig_c, esig_c,  q0_c,  q0_c);    f_ss_q0=scalar_contract(X_ss_q0, Fs);
        X_ss_q1=Xtensor(mpol, K, time, parity, esig_c, esig_c,  q1_c,  q1_c);    f_ss_q1=scalar_contract(X_ss_q1, Fs);
 


### 
        #f_ss_diff=f_ss_q0- f_ss_q1#Borrmann case 
        
        #For absorption anisotropy: Specify vector for polarization and beam direction (or pol and k for Borrmann)
        #(Ehkl, qhkl, rothkl)=(array([0,1,0]), array([-1,0,0]), array([1,1,0]) ); #hkl values for polarization, beam direction, psi rotation axis
        #(Ecart, qcart, rotcart)=(dot(B,Ehkl), dot(B,qhkl), dot(B,rothkl));  #convert to crystal Cartesian
        #next two lines not checked yet - can comment out these or whole absorption section
        #Urot=rotate_about_vector(rotcart, psi);
        #(Ecart, qcart)=(dot(Urot,Ecart), dot(Urot,qcart));  #rotate by psi about rothkl vector
        X_mu=Xtensor(mpol, K, time, parity, z_c, z_c, y_c , y_c);    f_0=scalar_contract(X_mu, Fs);#f_0 is absorption assuming that k is along y, pol is along z and sample is rotated about specified hkl vector (x) by azimuthal angle
        
        

    elif K==-1:
        #non-resonant magnetic scattering
        #spin and orbital componets (complex) for reflection are sk, lk 
        #BB and AA are B (spin) and A (orbit) coupling vectors from SWL, Blume etc
        #Need to check the signs...
        #Error in sign of third term. Old expression: BB=cross(e1,e0)+cross(q1_c,e1)*dot(q1_c,e0)+cross(q0_c,e0)*dot(q0_c,e1)-cross(cross(q1_c,e1),cross(q0_c,e0))
        
        e0=esig_c; e1=esig_c;
        BB=cross(e1,e0)+cross(q1_c,e1)*dot(q1_c,e0)-cross(q0_c,e0)*dot(q0_c,e1)-cross(cross(q1_c,e1),cross(q0_c,e0))
        AA=2*(1-dot(q0_c,q1_c))*cross(e1,e0)-cross(q0_c,e0)*dot(q0_c,e1)+cross(q1_c,e1)*dot(q1_c,e0)
        f_ss=1j*(dot(sk,BB)+dot(lk,AA)); 
 
        e0=esig_c; e1=epi1_c;
        BB=cross(e1,e0)+cross(q1_c,e1)*dot(q1_c,e0)-cross(q0_c,e0)*dot(q0_c,e1)-cross(cross(q1_c,e1),cross(q0_c,e0))
        AA=2*(1-dot(q0_c,q1_c))*cross(e1,e0)-cross(q0_c,e0)*dot(q0_c,e1)+cross(q1_c,e1)*dot(q1_c,e0)
        f_sp=1j*(dot(sk,BB)+dot(lk,AA));
        
        e0=epi0_c; e1=esig_c;
        BB=cross(e1,e0)+cross(q1_c,e1)*dot(q1_c,e0)-cross(q0_c,e0)*dot(q0_c,e1)-cross(cross(q1_c,e1),cross(q0_c,e0))
        AA=2*(1-dot(q0_c,q1_c))*cross(e1,e0)-cross(q0_c,e0)*dot(q0_c,e1)+cross(q1_c,e1)*dot(q1_c,e0)
        f_ps=1j*(dot(sk,BB)+dot(lk,AA));
        
        e0=epi0_c; e1=epi1_c;
        BB=cross(e1,e0)+cross(q1_c,e1)*dot(q1_c,e0)-cross(q0_c,e0)*dot(q0_c,e1)-cross(cross(q1_c,e1),cross(q0_c,e0))
        AA=2*(1-dot(q0_c,q1_c))*cross(e1,e0)-cross(q0_c,e0)*dot(q0_c,e1)+cross(q1_c,e1)*dot(q1_c,e0)
        f_pp=1j*(dot(sk,BB)+dot(lk,AA)); 

        f_0=0;   #for compatibility with tensor calc


    elif K==-2:
        #resonant magnetic scattering
        #Need to check the signs...
        e0=esig_c; e1=esig_c;
        E1E1=cross(e1,e0)
        f_ss=1j*(dot(e1k,E1E1));
 
        e0=esig_c; e1=epi1_c;
        E1E1=cross(e1,e0)
        f_sp=1j*(dot(e1k,E1E1));

        e0=epi0_c; e1=esig_c;
        E1E1=cross(e1,e0)
        f_ps=1j*(dot(e1k,E1E1));

        e0=epi0_c; e1=epi1_c;
        E1E1=cross(e1,e0)
        f_pp=1j*(dot(e1k,E1E1));

        f_0=0;   #for compatibility with tensor calc


    G=array([[f_ss,  f_ps], [f_sp,  f_pp]])

    #intensities for four main pol channels and totals
    Iss=abs(f_ss)**2; Isp=abs(f_sp)**2; Ips=abs(f_ps)**2; Ipp=abs(f_pp)**2; Is=Iss+Isp; Ip=Ips+Ipp;

    #Stokes parameters (incient beam) P3=lin, P1=45, P2=circ
    P3=1; P1=P2=0;###################+++++++++
    ##P3=0.0; P1=-0.1; P2=0.9
    mu=1./2*array([[1.+P3, P1-1.J*P2], [ P1+1.J*P2, 1.-P3]])
    
    eta=0*pi/4
    #Analyser matrix (reverse sign of eta compared to HoFe2 paper)
    A=array([[cos(eta), sin(eta)], [-cos(2*pol_theta)*sin(eta), cos(2*pol_theta)*cos(eta) ]])
    #intensity through polarizer
    I_0=dot(A, dot(G, dot(mu,dot(conjugate(G.T), conjugate(A.T))))).trace()
    
    eta=1*pi/4
    #Analyser matrix (reverse sign of eta compared to HoFe2 paper)
    A=array([[cos(eta), sin(eta)], [-cos(2*pol_theta)*sin(eta), cos(2*pol_theta)*cos(eta) ]])
    #intensity through polarizer
    I_45=dot(A, dot(G, dot(mu,dot(conjugate(G.T), conjugate(A.T))))).trace()
        
    eta=2*pi/4
    #Analyser matrix (reverse sign of eta compared to HoFe2 paper)
    A=array([[cos(eta), sin(eta)], [-cos(2*pol_theta)*sin(eta), cos(2*pol_theta)*cos(eta) ]])
    #intensity through polarizer
    I_90=dot(A, dot(G, dot(mu,dot(conjugate(G.T), conjugate(A.T))))).trace()
    
    eta=3*pi/4
    #Analyser matrix (reverse sign of eta compared to HoFe2 paper)
    A=array([[cos(eta), sin(eta)], [-cos(2*pol_theta)*sin(eta), cos(2*pol_theta)*cos(eta) ]])
    #intensity through polarizer
    I_135=dot(A, dot(G, dot(mu,dot(conjugate(G.T), conjugate(A.T))))).trace()
        
    #total intensity for given Stokes parameter
    I_tot=dot(G, dot(mu,conjugate(G.T))).trace()

    #circ pol total intensity
    #Stokes parameters (incient beam) P3=lin, P1=45, P2=circ
    P2=1; P1=P3=0;
    mu=1./2*array([[1.+P3, P1-1.J*P2], [ P1+1.J*P2, 1.-P3]])
    I_tot_circp=dot(G, dot(mu,conjugate(G.T))).trace()
    I_tot_circm=dot(G, dot(conjugate(mu),conjugate(G.T))).trace()


    #calculate intensity scaling for sample absorption
    Nhkl=array([0,0,1]);#surface normal hkl
    Nc=dot(B,Nhkl);
    q0c=dot(Uctheta,q0);
    q1c=dot(Uctheta,q1);
    sin_alpha=-dot(q0c,Nc);
    sin_beta=dot(q1c,Nc);
    #Iabsorb=dot(q1c,Nc)/(-dot(q0c,Nc)+dot(q1c,Nc));
    Iabsorb=(sin_alpha>0)*(sin_beta>0)*sin_beta/(sin_alpha+sin_beta);
    #mu_t=0.04;#best fit 303
    mu_t=0.01;
    atten_layer=exp(-mu_t*(1/sin_alpha+1/sin_beta)); #attenuation from surface layer
    Iabsorb=min(100,Iabsorb); #max enhancement 100
    Iabsorb=Iabsorb*atten_layer;    #correct for surface dead layer
    Iabslist+=[Iabsorb];
    Itotabs+=[I_tot*Iabsorb];#total intensity corrected for absorption

    psilist+=[psideg];  totlist+=[I_tot]; p0list+=[I_0]; p45list+=[I_45]; p90list+=[I_90]; p135list+=[I_135];circplist+=[I_tot_circp]; circmlist+=[I_tot_circm]; Isslist+=[Iss]; Isplist+=[Isp]; Ipslist+=[Ips]; Ipplist+=[Ipp]; Islist+=[Is]; Iplist+=[Ip]; f0_r+=[real(f_0)];f0_im+=[imag(f_0)];
    #print '%.0f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f' % (psideg, Iss, Isp, I_tot,Iabsorb,I_tot*Iabsorb)
    #print '%.0f\t%.4f\t%.4f\t%.4f' % (psideg,I_tot,Iabsorb,I_tot*Iabsorb)
    #print '%.0f\t%.4f\t%.4f\t%.4f\t%.4f' % (psideg, I_0, I_45, I_90, I_135 )
    #print '%.6f' % (I_tot*Iabsorb)
   
titlestr='hkl=[%.1f, %.1f, %.1f]   $\psi_0$=[%.1f, %.1f, %.1f]' % (tuple(hkl)+tuple(hkln)) 
   
#sig-sig, sig-pi, sig-total
plt.figure(1); plt.hold(True);
plt.plot(psilist, Islist, 'k',label='$\sigma$ Total',linewidth=2.0);
plt.plot(psilist, Isslist, 'r',label='$\sigma\sigma$',linewidth=2.0);
plt.plot(psilist, Isplist, 'b',label='$\sigma\pi$',linewidth=2.0); 
plt.hold(False);plt.legend(loc='best'); plt.ylabel('Intensity (aribtrary units)');plt.axis('tight'); plt.xlabel('$\psi$ (degrees)'); title(titlestr); grid(1)
savefig('/tmp/%s sigma.pdf' % filestring)


import pandas as pd
df=pd.DataFrame(zip(psilist,Isslist, Isplist, Islist))
df.to_csv('/media/sf_Dropbox/tmp/%s sigma.txt' % filestring, header=False, index=False, sep='\t')

#plt.show();

#pi-sig, pi-pi, pi-total
plt.figure(2); plt.hold(True);
plt.plot(psilist, Iplist, 'k',label='$\pi$ Total',linewidth=2.0);
plt.plot(psilist, Ipslist, 'r',label='$\pi\sigma$',linewidth=2.0);
plt.plot(psilist, Ipplist, 'b',label='$\pi\pi$',linewidth=2.0); 
plt.hold(False);plt.legend(loc='best'); plt.ylabel('Intensity (aribtrary units)');plt.axis('tight'); plt.xlabel('$\psi$ (degrees)'); title(titlestr)
#savefig('/media/sf_Dropbox/tmp/pi.pdf')
#plt.show();




#frac pi
fracpi=(array(Isplist)-array(Isslist))/(array(Isplist)+array(Isslist));
plt.figure(); plt.plot(psilist, fracpi, 'b',label='$(\pi-\sigma)/(\pi+\sigma)$',linewidth=2.0); plt.axis('tight'); plt.legend(loc='best'); plt.xlabel('$\psi$ (degrees)'); title(titlestr)
#plt.show();

#plt.figure(2); plt.plot(psilist, f0_r, 'r',psilist, f0_im,'g')
#plt.xlabel('psi (degrees)'); plt.ylabel('f_0 real(red) imag (green)');plt.axis('tight'); plt.show();



#plt.plot(psilist, Iabslist,'k',psilist,Itotabs,'r', psilist, totlist, 'b'  )    
#plt.plot(psilist,Itotabs,'r')    
plt.plot(psilist, p0list, 'b', psilist, p45list, 'g', p90list, 'r', psilist, p135list, 'm' ); plt.axis('tight');  plt.legend(('pol=0','pol=45','pol=90','pol=135'),loc='best'); plt.ylabel('Intensity (aribtrary units)'); plt.xlabel('$\psi$ (degrees)'); grid(1)
plt.show();
#plt.plot(psilist, totlist, 'r',psilist, pollist, 'b')
plt.plot(psilist, totlist, 'k-o'); plt.ylabel('Intensity (aribtrary units)')
#plt.hold(1)
#plt.plot(psilist, circplist, 'r',psilist, circmlist, 'b');
#plt.hold(0)
print '=== mean circp='+str(mean(circplist[0:-1]))
print '=== mean circm='+str(mean(circmlist[0:-1]))

#just type axis('tight') etc in iPython


print "=== Done"

#for ii in range(len(sglist)):
#    print ii,sglist[ii]


