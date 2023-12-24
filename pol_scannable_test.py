
#/dls_sw/i16/scripts/Users/steve/

#run('/dls_sw/i16/scripts/Users/steve/pol_scannable_test.py')

def trace2(array2):
    # trace of 2x2 array
    return (array2[0, 0] + array2[1, 1])

def dot(A, B):
    # dnp dot method fails with complex arrays
    Ar, Aim, Br, Bim = dnp.real(A), np.imag(A), np.real(B), np.imag(B)
    return dnp.dot(Ar, Br) - dnp.dot(Aim, Bim) + 1.J * dnp.dot(Ar, Bim) + 1.J * dnp.dot(Aim, Br) 

class polarizer(ScannableMotionBase):
    '''
    self.set_scatter_type - set scatter type; self.set_scatter_type('') to see options
    set_mag_vec - set magnetic vector to list/array
    returns:
        JA[0,0], JA[1,0], JA[0,1], JA[1,1]  - Polarization analyser Jones matrix components
        Idirect, Itot, Ipol                 - Direct beam intensity through analyser, total intensity from sample, polarization analyer intensity from sample
    '''

    def __init__(self,name):
        self.setName(name)
        self.setLevel(6)
        self.scatter_type = 'charge'    # default scattering mechanism
        self.setInputNames([])
        self.setExtraNames(['J00', 'J01', 'J10', 'J11', 'Idirect', 'Itot_' + self.scatter_type, 'Ipol_' + self.scatter_type])
        self.setOutputFormat(['%.4f','%.4f','%.4f','%.4f','%.3f', '%.3f', '%.3f'])
        self.m1, self.m2, self.m3 = 0., 0., 1.

        
    def set_scatter_type(self, scatter_type):
        scatter_types = ['charge', 'magE1E1', 'magSpin'] # scattering mechanism
        if scatter_type in scatter_types:
          self.scatter_type = scatter_types
        else:
          print('=== scatter type must be one of: ', scatter_types)
        print('Scatter type is: ', self.scatter_type)
        
    def set_mag_vec(self, magvec):
        [self.m1, self.m2, self.m3] =  magvec

    def BB(q, q1, e, e1)  # spin scattering vector
        dnp.cross(e1, e) + dnp.cross(q1, e1)*dnp.dot(q1, e) 
        - dnp.cross(q, e)*dnp.dot(q, e1) - dnp.cross(dnp.cross(q1, e1),dnp.cross(q, e))
       

    def getPosition(self):
        # get motor positions
        xi_deg = stokes()
        tthp_deg = 2 * thp() #assume thp is correct; tthp may be offset
        gamma_deg = gam()
        delta_deg = delta()
            
        P1, P2, P3 = 1., 0., 0.
        _mu = dnp.array([[P1/2. + 1/2., P2/2. - 1J*P3/2.], [P2/2. + 1J*P3/2., 1/2. - P1/2.]])
    
            
        xi_rad, tthp_rad, gamma_rad, delta_rad = xi_deg * dnp.pi/180, tthp_deg * dnp.pi/180, gamma_deg * dnp.pi/180, delta_deg * dnp.pi/180
            
        # calculate J and Icharge
            
        self.JA = dnp.array([[dnp.cos(xi_rad), -dnp.sin(xi_rad)], [dnp.sin(xi_rad)*np.cos(tthp_rad), np.cos(tthp_rad)*np.cos(xi_rad)]])      

        Rgammadelta = dnp.array([[dnp.cos(gamma_rad), -dnp.sin(delta_rad)*dnp.sin(gamma_rad), dnp.sin(gamma_rad)*dnp.cos(delta_rad)], 
                                [0., dnp.cos(delta_rad), dnp.sin(delta_rad)], 
                                [-dnp.sin(gamma_rad), -dnp.sin(delta_rad)*dnp.cos(gamma_rad), dnp.cos(delta_rad)*dnp.cos(gamma_rad)]])

        x, y, z = dnp.array([1., 0., 0.]), dnp.array([0., 1., 0.]), dnp.array([0., 0., 1.])
        x1, y1, z1 = dnp.dot(Rgammadelta,x), dnp.dot(Rgammadelta,y), dnp.dot(Rgammadelta,z)
        
        U3 = (z - z1) / dnp.linalg.norm(z - z1)
        U1 = (z + z1) / dnp.linalg.norm(z + z1)
        U2 = dnp.cross(U3, U1)
        
        m = self.m1 * U1 + self.m2 * U2 + self.m3 * U3

        if self.scatter_type == 'charge':
            self.J = dnp.array([[dnp.cos(gamma_rad), 0], [-dnp.sin(delta_rad)*dnp.sin(gamma_rad), dnp.cos(delta_rad)]])
        elif self.scatter_type == 'magE1E1':
            self.J =  dnp.array([
              [-m[1] * dnp.sin(gamma_rad), m[0] * dnp.sin(gamma_rad) + m[2] * dnp.cos(gamma_rad)],
              [-m[1] * dnp.sin(delta_rad)*dnp.cos(gamma_rad) - m[2] * dnp.cos(delta_rad), m[0] * dnp.sin(delta_rad)*dnp.cos(gamma_rad) - m[2] * dnp.sin(delta_rad)*dnp.sin(gamma_rad)]])
        elif self.scatter_type == 'magSpin':
            self.J =  dnp.array([
                [dnp.dot(m, BB(q, q1, x, x1)), dnp.dot(m, BB(q, q1, y, x1))],
                [dnp.dot(m, BB(q, q1, x, y1)), dnp.dot(m, BB(q, q1, y, y1))]
            ])
        else:
            raise('Scatter type', self.scatter_type, ' not supported')

            

        
        self.Itot = trace2(dot(self.J, dot(_mu, dot(dnp.conjugate(self.J).T))))
        #result should be real but convert to allow scan
        self.Itot =  float(dnp.real(self.Idirect))


        self.Ipol = trace2(dot(self.JA, dot(self.J, dot(_mu, dot(dnp.conjugate(self.J).T, dnp.conjugate(self.JA).T)))))
        #result should be real but convert to allow scan
        self.Ipol =  float(dnp.real(self.Idirect))
        
        self.Idirect = trace2(dot(self.JA, dot(_mu, dnp.conjugate(self.JA).T)))
        #result should be real but convert to allow scan
        self.Idirect =  float(dnp.real(self.Idirect))
        
        
        return [self.JA[0,0], self.JA[1,0], self.JA[0,1], self.JA[1,1], self.Idirect, self.Itot, self.Ipol]

    def isBusy(self):
        return 0     

#error in this:
#dnp.trace(dnp.matmul(ss.JA, dnp.matmul(ss.JC, dnp.matmul(mu, dnp.matmul(dnp.conjugate(ss.JC).T, dnp.conjugate(ss.JA).T)))))
        
#### test out dnp maths functions

pp = polarizer('I16 PA')



pp()
pp


# dnp.dot(self.JA, 
#         dnp.dot(self.JC, 
#                 dnp.dot(_mu, 
#                         dnp.dot(dnp.conjugate(self.JC).T, 
#                                 dnp.conjugate(self.JA).T))))

# no conjugate - real

P1, P2, P3 = 1., 0., 0.
_mu = dnp.array([[P1/2. + 1/2., P2/2. - 1J*P3/2.], [P2/2. + 1J*P3/2., 1/2. - P1/2.]])

print(); print(_mu); print()
jc_mu_jc = dnp.dot(pp.JC, dnp.dot(_mu, pp.JC.T))
print(jc_mu_jc); print()
ja_jc_mu_jc_ja = dnp.dot(pp.JA, dnp.dot(jc_mu_jc, pp.JA.T))
print(ja_jc_mu_jc_ja); print()

mm1 = dnp.array([[1., 0.], [0., 0.]])
mm2 = dnp.array([[0.70712579, 0.70708777],
 [-0.70708778, 0.70712578]])
mm3 = dnp.dot(mm1, mm2)

# debug on I16
# add charge scattering
# add magnetic scattering


#19/6/23 tests

#P1, P2, P3 = 1, 0, 0
#mu = dnp.array([[P1/2. + 1/2., P2/2. - 1J*P3/2.], [P2/2. + 1J*P3/2., 1/2 - P1/2.]]) #ok
#dnp.dot(mu, mu) #ok
#dnp.trace(mu) #fail: trace(): expected 1 or 4 args; got 2
#mu.T #ok
#dnp.conjugate(mu)    #ok


#>>> ss()
#[0.9999999999897039, 4.5378560551696824e-06, -4.5378560551696824e-06, 0.9999999999897039, 0j]
#>>> ss
#I16 PA : J00: 1.0000 J01: 0.0000 J10: -0.0000 J11: 1.0000 Icharge: 0j

#18/7/23

#>>> tthp
#tthp : -112.13deg mot(-100.20:209.80)
#>>> 
#>>> gam
#gam : 0.00165 (-1.00000:120.00000)
#>>> delta
#delta : 28.90050 (-1.35000:109.65000)


#>>> pp
#I16 PA : J00: 1.0000 J01: -0.0000 J10: 0.0000 J11: 1.0000 Icharge: (0.999999999143-1.39197306961e-05j)
#>>> pos stokes 45
#Move completed: stokes : 45.000
#>>> pp
#I16 PA : J00: 0.7071 J01: 0.7071 J10: -0.7071 J11: 0.7071 Icharge: (0.50000870478-6.95998652219e-06j)
#>>> pos stokes 90
#Move completed: stokes : 90.000
#>>> pp
#I16 PA : J00: 0.0000 J01: 1.0000 J10: -1.0000 J11: 0.0000 Icharge: (4.03552453533e-10-5.60360069092e-15j)

# >>> _mu
# array([[1.0000000 + 0.0000000j, 0.0000000 + 0.0000000j],
#  [0.0000000 + 0.0000000j, 0.0000000 + 0.0000000j]])
# >>> pp.JA
# array([[0.70712579, -0.70708778],
#  [0.70708777, 0.70712578]])
# >>> dnp.dot(_mu, pp.JA)
# array([[0.70712579 - 0.70708778j, 0.0000000 + 0.0000000j],
#  [0.0000000 + 0.0000000j, 0.0000000 + 0.0000000j]])


# error in dnp (_mm3 shold be real)
_mm1 = dnp.array([[1.0000000 + 0.0000000j, 0.0000000 + 0.0000000j],
 [0.0000000 + 0.0000000j, 0.0000000 + 0.0000000j]])
_mm2 = dnp.array([[0.70712579, -0.70708778],
 [0.70708777, 0.70712578]])
_mm3 = dnp.dot(_mm1, _mm2)
#array([[0.70712579 - 0.70708778j, 0.0000000 + 0.0000000j],
# [0.0000000 + 0.0000000j, 0.0000000 + 0.0000000j]])

'''

    def NonResonantMagneticScatteringMatrix(self, sk, lk, esig_c, e0pi_c, e1pi_c, q0_c,  q1_c):
        '''
        Calculate 2x2 scattering amplitude matrix for non-resonant magnetic scattering
        spin and orbital components (complex) for reflection are sk, lk 
        BB and AA are B (spin) and A (orbit) coupling vectors from SWL, Blume etc
        '''
        e0=esig_c; e1=esig_c;
        BB=np.cross(e1,e0)+np.cross(q1_c,e1)*np.dot(q1_c,e0)-np.cross(q0_c,e0)*np.dot(q0_c,e1)-np.cross(np.cross(q1_c,e1),np.cross(q0_c,e0))
        AA=2*(1-np.dot(q0_c,q1_c))*np.cross(e1,e0)-np.cross(q0_c,e0)*np.dot(q0_c,e1)+np.cross(q1_c,e1)*np.dot(q1_c,e0)
        f_ss=1j*(np.dot(sk,BB)+np.dot(lk,AA)); 
 
        e0=esig_c; e1=e1pi_c;
        BB=np.cross(e1,e0)+np.cross(q1_c,e1)*np.dot(q1_c,e0)-np.cross(q0_c,e0)*np.dot(q0_c,e1)-np.cross(np.cross(q1_c,e1),np.cross(q0_c,e0))
        AA=2*(1-np.dot(q0_c,q1_c))*np.cross(e1,e0)-np.cross(q0_c,e0)*np.dot(q0_c,e1)+np.cross(q1_c,e1)*np.dot(q1_c,e0)
        f_sp=1j*(np.dot(sk,BB)+np.dot(lk,AA));
        
        e0=e0pi_c; e1=esig_c;
        BB=np.cross(e1,e0)+np.cross(q1_c,e1)*np.dot(q1_c,e0)-np.cross(q0_c,e0)*np.dot(q0_c,e1)-np.cross(np.cross(q1_c,e1),np.cross(q0_c,e0))
        AA=2*(1-np.dot(q0_c,q1_c))*np.cross(e1,e0)-np.cross(q0_c,e0)*np.dot(q0_c,e1)+np.cross(q1_c,e1)*np.dot(q1_c,e0)
        f_ps=1j*(np.dot(sk,BB)+np.dot(lk,AA));
        
        e0=e0pi_c; e1=e1pi_c;
        BB=np.cross(e1,e0)+np.cross(q1_c,e1)*np.dot(q1_c,e0)-np.cross(q0_c,e0)*np.dot(q0_c,e1)-np.cross(np.cross(q1_c,e1),np.cross(q0_c,e0))
        AA=2*(1-np.dot(q0_c,q1_c))*np.cross(e1,e0)-np.cross(q0_c,e0)*np.dot(q0_c,e1)+np.cross(q1_c,e1)*np.dot(q1_c,e0)
        f_pp=1j*(np.dot(sk,BB)+np.dot(lk,AA)); 
        
        
        self.G=np.array([[f_ss,  f_ps], [f_sp,  f_pp]])    #scattering matrix
        return self.G



    def E1E1ResonantMagneticScatteringMatrix(self, mk, esig_c, e0pi_c, e1pi_c, q0_c,  q1_c):
        '''
        Calculate 2x2 scattering amplitude matrix for E1E1 resonant magnetic scattering
        '''
        e0=esig_c; e1=esig_c;
        E1E1=np.cross(e1,e0)
        f_ss=1j*(np.dot(mk,E1E1));
 
        e0=esig_c; e1=e1pi_c;
        E1E1=np.cross(e1,e0)
        f_sp=1j*(np.dot(mk,E1E1));

        e0=e0pi_c; e1=esig_c;
        E1E1=np.cross(e1,e0)
        f_ps=1j*(np.dot(mk,E1E1));

        e0=e0pi_c; e1=e1pi_c;
        E1E1=np.cross(e1,e0)
        f_pp=1j*(np.dot(mk,E1E1));
        
        self.G=np.array([[f_ss,  f_ps], [f_sp,  f_pp]])    #scattering matrix
        return self.G




