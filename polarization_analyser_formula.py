import numpy as np
import quickfit

#I_str='np.real( scalefac*(((StokesP1/2. - 0.5*1j*StokesP2)*np.cos(x*np.pi/180) + (-StokesP3/2. + 1/2.)*np.sin(x*np.pi/180))*np.sin(x*np.pi/180) + ((StokesP1/2. + 0.5*1j*StokesP2)*np.sin(x*np.pi/180) + (StokesP3/2. + 1/2.)*np.cos(x*np.pi/180))*np.cos(x*np.pi/180) + (-(StokesP1/2. - 0.5*1j*StokesP2)*np.sin(x*np.pi/180)*np.cos(2*theta_pol*np.pi/180) + (-StokesP3/2. + 1/2.)*np.cos(2*theta_pol*np.pi/180)*np.cos(x*np.pi/180))*np.cos(2*theta_pol*np.pi/180)*np.cos(x*np.pi/180) - ((StokesP1/2. + 0.5*1j*StokesP2)*np.cos(2*theta_pol*np.pi/180)*np.cos(x*np.pi/180) - (StokesP3/2. + 1/2.)*np.sin(x*np.pi/180)*np.cos(2*theta_pol*np.pi/180))*np.sin(x*np.pi/180)*np.cos(2*theta_pol*np.pi/180)) )'
I_str='scalefac*((StokesP1*np.sin(x*np.pi/180)/2. + (StokesP3/2. + 1/2.)*np.cos(x*np.pi/180))*np.cos(x*np.pi/180) + (StokesP1*np.cos(x*np.pi/180)/2. + (-StokesP3/2. + 1/2.)*np.sin(x*np.pi/180))*np.sin(x*np.pi/180) + (-StokesP1*np.sin(x*np.pi/180)*np.cos(2*theta_pol*np.pi/180)/2. + (-StokesP3/2. + 1/2.)*np.cos(2*theta_pol*np.pi/180)*np.cos(x*np.pi/180))*np.cos(2*theta_pol*np.pi/180)*np.cos(x*np.pi/180) - (StokesP1*np.cos(2*theta_pol*np.pi/180)*np.cos(x*np.pi/180)/2. - (StokesP3/2. + 1/2.)*np.sin(x*np.pi/180)*np.cos(2*theta_pol*np.pi/180))*np.sin(x*np.pi/180)*np.cos(2*theta_pol*np.pi/180))'


PA=quickfit.fit_func('PA pol/eta scan intensity (angles in degrees)', ['theta_pol', 'scalefac', 'StokesP1', 'StokesP3'], I_str)
PA.params['theta_pol'].vary=False




'''
from sympy import *
#Stokes parameters (incient beam) P3=lin, P1=45, P2=circ
(P1,P2,P3, scalefac)=symbols('StokesP1,StokesP2,StokesP3, scalefac',real=True)
(theta_pol,x)=symbols('theta_pol,x',real=True) #x is eta_pol
P2=0 # measure only linear pol P1, P3; P2 calculated as being 'all that is not linear'
mu=Matrix([[1+P3, P1-1J*P2], [ P1+1J*P2, 1-P3]])/2
#Analyser matrix
A=Matrix([[cos(x), sin(x)], [-cos(2*theta_pol)*sin(x), cos(2*theta_pol)*cos(x) ]])
I_pol=scalefac*( A*mu*A.T.conjugate()).trace()            #intensity through polarizer
I_pol.simplify()    #formula for fit
I_str=str(I_pol)

#I_str=I_str.replace('I', '1j') #formula for fit
I_str=I_str.replace('sin','np.sin')
I_str=I_str.replace('cos','np.cos')
I_str=I_str.replace('/2','/2.')
I_str=I_str.replace('theta_pol','theta_pol*np.pi/180')
I_str=I_str.replace('x','x*np.pi/180')
#print "I_str='np.real(", I_str, ")'" #paste this line above and comment out sympy section
print "I_str='"+I_str+"'" #paste this line above and comment out sympy section
'''