#23/4/14 save old versions of tools & examples before attempting to incorporate magnetic symmetry

import sys, pprint
from numpy import *
from numpy.linalg import *
from numpy.random import *
from scipy.misc import *
sys.path.append('/home/spc93/python/PyCifRW-3.1.2')
sys.path.append('/media/DCF0769CF0767D18/python/PyCifRW-3.1.2')
import CifFile
import matplotlib.pyplot as plt
set_printoptions(precision=3, suppress=True)


#from interactive import diag
from copy import copy, deepcopy

def genpos2matvec(gen_pos_string):
    'convert general position string to vector/matrix form (floats) using lists as row vectors'
    #gp=gen_pos_string
    gp=gen_pos_string.lower(); # 27/9/16 modified to convert to lower for CIF files that use capital XYZ
    x=y=z=0.; vec=list(eval(gp.replace('/','./')))
    x=1.; y=z=0.; m0=list(eval(gp.replace('/','./'))); m0[0]=m0[0]-vec[0]; m0[1]=m0[1]-vec[1];m0[2]=m0[2]-vec[2];
    y=1.; x=z=0.; m1=list(eval(gp.replace('/','./'))); m1[0]=m1[0]-vec[0]; m1[1]=m1[1]-vec[1];m1[2]=m1[2]-vec[2];
    z=1.; x=y=0.; m2=list(eval(gp.replace('/','./'))); m2[0]=m2[0]-vec[0]; m2[1]=m2[1]-vec[1];m2[2]=m2[2]-vec[2];
    return [array([m0, m1, m2]).T, array(vec)]

def invert_spacegroup(sg):
    #perform inversion on spacegroup operators by reversing vector parts
    #seems to change original group as well as copy - OK with deepcopy
    newsg=deepcopy(sg)
    for sgop in newsg:
        sgop[1]=-sgop[1]
    return newsg

def T_symmetry_consistent_with_fm_vector(sg,magvec,Bmat=array([[1, 0, 0],[0,1,0],[0,0,1]])):
    #sg must already contain third field (time symmetry) although it is not used
    #apply matrix part of symmetry operator
    #magvec in crystal cartesian coordinates - need to convert first if not
    #tries to adjust magnetic symmetry so vectors are exactly consistent with symmetry
    errmax=1.e-6
    print '~~~~~~~~~~'; print 'magvec'; print magvec; print '~~~~~~~~~~'; print 'transformed vectors'; ######################delete line
    Tsym=[]
    for symop in sg:
        magvecnew=transform_cart(magvec, crystal_to_cart_operator(symop[0], Bmat),P=+1)
        print magvecnew######################delete line
        if sqrt((magvecnew-magvec).dot(magvecnew-magvec))<errmax:
            Tsym+=[1]  #magvec is symmetric under transformation
        elif sqrt((magvecnew+magvec).dot(magvecnew+magvec))<errmax:
            Tsym+=[-1]  #magvec is antisymmetric under transformation   
        else:
            Tsym+=[0]  #magvec neither symmetric nor antisymmetric   
    return Tsym

def T_symmetry_favouring_fm_vector(sg,magvec,Bmat=array([[1, 0, 0],[0,1,0],[0,0,1]])):
    #sg must already contain third field (time symmetry) although it is not used
    #apply matrix part of symmetry operator
    #magvec in crystal cartesian coordinates - need to convert first if not
    #changes time sym if necessary to make mag vector more parallel than antiparallel
    errmax=1.e-6
    print '~~~~~~~~~~'; print 'magvec'; print magvec; print '~~~~~~~~~~'; print 'transformed vectors'; ######################delete line
    Tsym=[]
    for symop in sg:
        magvecnew=transform_cart(magvec, crystal_to_cart_operator(symop[0], Bmat),P=+1)
        print 'new mag vec', magvecnew######################delete line
        if ((magvecnew).dot(magvec))>errmax:
            Tsym+=[1]
        elif ((magvecnew).dot(magvec))<-errmax:
            Tsym+=[-1]#reverse vector
        else:
            Tsym+=[0]# neither + nor -
    return Tsym




def firstCell(V):
    #fold V back to first unit cell (0..1)
    #return array([z-int(z)+int(z<0) for z in V])
    return array([z-floor(z) for z in V])

def isGroup(G):
    tol=0.0000001
    #group is a list of [mat, vec, timescalar]
    eye_index=-1
    for ind in range(len(G)):
        if all(abs(G[ind][0]-eye(3))<tol) and all(abs(G[ind][1]-zeros(3))<tol) and abs(G[ind][2]-1)<tol:
            eye_index=ind
    if eye_index !=0:
        print '=== Warning: Identity not first element'
    for S1 in G:
        for S2 in G:
            M3=dot(S1[0], S2[0])
            #V3=firstCell(S1[1] + S2[1]);    #fold back to first unit cell
            V3=firstCell(S1[1] + dot(S1[0],S2[1]));    #fold back to first unit cell
            T3=S1[2] * S2[2]
            n=0
            for S3 in G:
                if all(abs(M3-S3[0])<tol) and all(abs(V3-firstCell(S3[1]))<tol) and abs(T3-S3[2])<tol:
                    n+=1
            if n!=1:
                print '=== Warning: Not a group!'
                print '=== There should be one occurence of the following symmetry operator but were %i' % n
                print M3, V3, T3
                print '=== Derived from'
                print S1
                print '=== and'
                print S2
                return False
    return True

def set_spacegroup_timesym(sg, timesym):
    for ii in range(len(sg)):
        sg[ii][2]=timesym[ii]
    return sg
    
# def mag_group_consistent_with_fm_vector(sg,magvec,Bmat=array([[1, 0, 0],[0,1,0],[0,0,1]])):
#     #sg must already contain thirst field (time symmetry) although it is not used
#     #apply matrix part of symmetry operator
#     #magvec in crystal cartesian coordinates - need to convert first if not
#     errmax=1.e-6
#     newgroup=copy(sg)
#     print '~~~~~~~~~~'; print 'magvec'; print magvec; print '~~~~~~~~~~'; print 'transformed vectors'; ######################delete line
#     for symop in newgroup:
#         magvecnew=transform_cart(magvec, crystal_to_cart_operator(symop[0], Bmat),P=+1)
#         print magvecnew######################delete line
#         if sqrt((magvecnew-magvec).dot(magvecnew-magvec))<errmax:
#             symop[2]=1  #magvec is symmetric under transformation
#         elif sqrt((magvecnew+magvec).dot(magvecnew+magvec))<errmax:
#             symop[2]=-1  #magvec is antisymmetric under transformation   
#         else:
#             symop[2]=0  #magvec neither symmetric nor antisymmetric   
#     print '~~~~~~~~~~';######################delete line
#     return newgroup
#     
    
    
#mag = np.sqrt(x.dot(x))    
# def apply_sym(T, symop_list, Bmat=array([[1, 0, 0],[0,1,0],[0,0,1]]), P=0):
#     #apply point sym ops in symop_list to tensot T of rank K
#     #Optional Bmat is used to transform arrays to Cartesian from crystal basis
#     Tnew=T*0.0
#     for sym in symop_list:
#         Tnew+=transform_cart(T, crystal_to_cart_operator(sym, Bmat),P)
#     return Tnew
    
# def spacegroup_list_from_genpos_list(genposlist):
#     sglist=[];
#     for genpos in genposlist:
#         #sglist+=[genpos2matvec(genpos)]
#         sglist+=[genpos2matvec(genpos)+[1]] #add +1 to indicate time symmetry
#     return sglist    
    

def crystal_to_cart_operator(S, B):
    pass
    #transform crystal sym op to Cart sym op using B matrix
    Snew=dot(inv(B.T),  dot(S, B.T));
    return Snew

def rotate_about_vector(vec, ang):
    tol=0.999999
    #Generate unitary matrix that rotates about vector vec by angle ang (rad)
    #Uvz rotates vec onto z (assume it leaves (vec x z) unchanged)
    #check first that the rotation axis is not z (else just return z rotation matrix to avoid divide by zero)
    z=array([0,0,1])
    vec=array(vec)/norm(array(vec)); # in case its not already a unit vector array
    if dot(z,vec)>tol:  #parallel
        return array([[cos(ang), -sin(ang), 0],[sin(ang), cos(ang), 0],[0, 0, 1]]);
    elif dot(z,vec)<-tol:    #antiparallel
        return array([[cos(ang), sin(ang), 0],[-sin(ang), cos(ang), 0],[0, 0, 1]]);
    else:
        m1=array([vec,cross(vec,z)/norm(cross(vec,z)),cross(vec,cross(vec,z))/norm(cross(vec,cross(vec,z)))]).transpose()  #initial matrix (set of vectors); transpose as the vectors are row vectors in np and need to be collumn vectors
        m2=array([z,cross(vec,z)/norm(cross(vec,z)),cross(z,cross(vec,z))/norm(cross(z,cross(vec,z)))]).transpose()      #rotated matrix (set of vectors)
        Uvz=dot(m2, inv(m1))
        Rz=array([[cos(ang), -sin(ang), 0],[sin(ang), cos(ang), 0],[0, 0, 1]]);  #rotate about z
        U=dot(inv(Uvz),dot(Rz, Uvz))
        return U

def gen_twelve_fold_about_v(v):
    #symmetry group for 12-fold rotation about v (used for external vector/axial field)
    m=rotate_about_vector(v, pi/6); 
    group=[eye(3)];
    for i in range(11):
        last=group[-1]
        group+=[dot(last,m)]
    return group

def intersection(spacegroup, pointgroup):
    #calculate intersection of spacegroup with a pointgroup (e.g. an external field)
    #returns subset of spacegroup that has matrices in pointgroup
    newgroup=[]
    for sg in spacegroup:
        for pg in pointgroup:
            if allclose(sg[0], pg, atol=1e-6):
                newgroup+=[sg]
                break
    return newgroup

def transform_cart_1(T, S):
    #transform Cart tensor rank K using symmetry operator S
    #doesn't seem to work - come back to this (doesn't preserve trace if S is unitary)
    k=len(T.shape); #rank of T
    tnew=T*1.0;     #create copy of T
    for kk in range(k):
        tnew=tensordot(tnew, S, axes=(kk, -1))    #contract each index in turn
    return tnew


def transform_cart(T, S, P=0):
    #transform Cart tensor rank K using symmetry operator S
    #If optional parameter P (parity) is given then a correction is made to account for the otherwise incorrect
    #tranformation of Cartesian tensors derived from spherical pseudotensors (see Mittelwihr paper)
    d=3
    k=len(T.shape); #rank of T
    if P==0:
        Sfac=1
    elif P==1 or P==-1:
        Sfac=det(S)**((3+P*(-1)**k)/2)
    else:
        raise ValueError('Parity should be +1 (even), -1 (odd) or 0 (ignored)')
        
    ##### delete next two lines - diagnostics only
    if Sfac==-1:
        print "===Applying sign change for pseudotensor transormation"
    
    tnew=T*0.0;
    if k==0:
        tnew=T
    elif k==1:
        for i in range(d):
            for ii in range(d):
                tnew[i]+=S[i, ii]*T[ii]
    elif k==2:
        for i in range(d):
            for ii in range(d):
                for j in range(d):
                    for jj in range(d):
                        tnew[i, j]+=S[i, ii]*S[j, jj]*T[ii, jj]
    elif k==3:
        for i in range(d):
            for ii in range(d):
                for j in range(d):
                    for jj in range(d):
                        for k in range(d):
                            for kk in range(d):
                                tnew[i, j, k]+=S[i, ii]*S[j, jj]*S[k, kk]*T[ii, jj, kk]
    elif k==4:
        for i in range(d):
            for ii in range(d):
                for j in range(d):
                    for jj in range(d):
                        for k in range(d):
                            for kk in range(d):
                                for l in range(d):
                                    for ll in range(d):
                                        tnew[i, j, k, l]+=S[i, ii]*S[j, jj]*S[k, kk]*S[l, ll]*T[ii, jj, kk, ll]                             
    elif k==5:
        for i in range(d):
            for ii in range(d):
                for j in range(d):
                    for jj in range(d):
                        for k in range(d):
                            for kk in range(d):
                                for l in range(d):
                                    for ll in range(d):
                                        for m in range(d):
                                             for mm in range(d): 
                                                 tnew[i, j, k, l, m]+=S[i, ii]*S[j, jj]*S[k, kk]*S[l, ll]*S[m, mm]*T[ii, jj, kk, ll, mm]                                        
    elif k==6:
        for i in range(d):
            for ii in range(d):
                for j in range(d):
                    for jj in range(d):
                        for k in range(d):
                            for kk in range(d):
                                for l in range(d):
                                    for ll in range(d):
                                        for m in range(d):
                                             for mm in range(d): 
                                                 for n in range(d):
                                                     for nn in range(d): 
                                                         tnew[i, j, k, l, m, n]+=S[i, ii]*S[j, jj]*S[k, kk]*S[l, ll]*S[m, mm]*S[n, nn]*T[ii, jj, kk, ll, mm, nn]                                        
    else:
        raise ValueError('Tranformation for this tensor rank not coded: K=',k)
    return tnew*Sfac   #Apply correction factor Sfac (+/-1)


def apply_sym(Tensor, symop_list, Bmat=array([[1, 0, 0],[0,1,0],[0,0,1]]), P=None, T=+1):
    #apply point sym ops in symop_list to tensor of rank K
    #Optional Bmat is used to transform arrays to Cartesian from crystal basis
    #Default time (T) sym +1; no default for parity (P)
    Tnew=Tensor*0.0
    for sym in symop_list:
        tsign=T**((1-sym[1])/2) #sign change of time-odd tensor under time inversion
        print '=== in apply_sym: T, sym[1], signchange=',T, sym[1], tsign ########### delete #############
        Tnew+=transform_cart(Tensor, crystal_to_cart_operator(sym[0], Bmat),P)*tsign
    return Tnew


def apply_sym_old(Tensor, symop_list, Bmat=array([[1, 0, 0],[0,1,0],[0,0,1]]), P=0):
    #apply point sym ops in symop_list to tensot T of rank K
    #Optional Bmat is used to transform arrays to Cartesian from crystal basis
    Tnew=Tensor*0.0
    for sym in symop_list:
        Tnew+=transform_cart(Tensor, crystal_to_cart_operator(sym, Bmat),P)
    return Tnew





def calc_SF(Tensor, R, hkl, spacegroup_list, Bmat=array([[1, 0, 0],[0,1,0],[0,0,1]]), P=None, T=+1):
    #calc structure factor tensor for symop_list to tensot T of rank K at position R hkl=Q
    #Optional Bmat is used to transform arrays to Cartesian from crystal basis
    Tnew=Tensor*0.0
    for sym in spacegroup_list:
        mat=sym[0]
        vec=sym[1]
        time=sym[2]
        tsign=T**((1-time)/2) #sign change of time-odd tensor under time inversion
        newR=dot(mat, R)+vec
        phase=exp(pi*2.j * dot(hkl, newR))
        #print 'phase:',phase
        #print transform_cart(T, mat)
        #Tnew+=transform_cart(T, mat)*phase
        newbit=transform_cart(Tensor, crystal_to_cart_operator(mat, Bmat),P)*phase*tsign
        #print 'symop (crys, cart):', mat, crystal_to_cart_operator(mat, Bmat)
        #print 'crytocart:',crystal_to_cart_operator(mat, Bmat)[0,0]
        #print 'transcart:',transform_cart(T, crystal_to_cart_operator(mat, Bmat), P)[0,0]
        #print 'newbit:',newbit[0,0]
        Tnew=Tnew+newbit
        #print 'Tnew:',Tnew[0,0]
        #raise ValueError('for debug')
    return Tnew

def calc_scalar_SF_from_sites(sites, hkl):
    #calculate scalar structure factor from set of equivalent atom coordinates
    F=0
    for site in sites:
        F+=exp(pi*2.j * dot(hkl, site))
    return F




def calc_SF_old(T, R, hkl, spacegroup_list, Bmat=array([[1, 0, 0],[0,1,0],[0,0,1]]), P=0):
    #calc structure factor tensor for symop_list to tensot T of rank K at position R hkl=Q
    #Optional Bmat is used to transform arrays to Cartesian from crystal basis
    Tnew=T*0.0
    for sym in spacegroup_list:
        mat=sym[0]
        vec=sym[1]
        newR=dot(mat, R)+vec
        phase=exp(pi*2.j * dot(hkl, newR))
        #print 'phase:',phase
        #print transform_cart(T, mat)
        #Tnew+=transform_cart(T, mat)*phase
        newbit=transform_cart(T, crystal_to_cart_operator(mat, Bmat),P)*phase
        #print 'symop (crys, cart):', mat, crystal_to_cart_operator(mat, Bmat)
        #print 'crytocart:',crystal_to_cart_operator(mat, Bmat)[0,0]
        #print 'transcart:',transform_cart(T, crystal_to_cart_operator(mat, Bmat), P)[0,0]
        #print 'newbit:',newbit[0,0]
        Tnew=Tnew+newbit
        #print 'Tnew:',Tnew[0,0]
        #raise ValueError('for debug')
    return Tnew


#def SF_symmetry(R, hkl, spacegroup_list, Bmat=array([[1, 0, 0],[0,1,0],[0,0,1]])):
#    #Look at symmetry of SF
#    #Optional Bmat is used to transform arrays to Cartesian from crystal basis
#    Podd=0;     #assumed not parity odd
#    Peven=0;    #assumed not parity even
#    Allowed=1;  #allowed untill proved otherwise
#    inv_mat=array([[-1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,-1.0]]); #inversion operator
#    identity_mat=array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]); #identity operator
#    identy_phase=exp(pi*2.j * dot(hkl, R)); #phase for initial site vector
#    #print 'identiy phase', identy_phase
#    for sym in spacegroup_list:
#
#        mat=sym[0]
#        vec=sym[1]
#        newR=dot(mat, R)+vec
#        phase=exp(pi*2.j * dot(hkl, newR))/identy_phase;    #change phases so first one is unity
#        #print mat, vec, phase
#
#        #delete next two lines
#        #if allclose(mat,identity_mat, rtol=.001):
#        #    print 'identity operator has phase '+ str(phase)
#
#        if allclose(mat,inv_mat, rtol=.001):
#            #print 'inversion operator has phase '+ str(phase)
#            if allclose(phase,-1.0):
#                Podd=1;     #it is parity odd
#            if allclose(phase,1.0):
#                Peven=1;    #it is parity even
#
#        if allclose(mat,identity_mat, rtol=.001):
#            #print 'inversion operator has phase '+ str(phase)
#            if allclose(phase,-1.0):
#                Allowed=0;     #all tensors will vanish
#    if Podd==1 and Peven==1:
#        Allowed=0
#    if Allowed==0:
#        print '=== No scattering tensor can exist'
#    else:
#        if Podd==1 and Peven==0:
#            print '=== Scattering tensor is parity odd'
#        elif Podd==0 and Peven==1:
#            print '=== Scattering tensor is parity even'
#        elif Podd==0 and Peven==0:
#            print '=== Scattering tensor is of mixed parity'
#    return [Allowed, Peven, Podd]


def SF_symmetry(R, hkl, spacegroup_list):
    #analyse symmetry of any possible structure factor (mainly for information)
    #returns [sym_phases, gen_scalar_allowed, site_scalar_allowed, tensor_allowed, Psym, Tsym, PTsym]
    tol=1e-6
    inv_mat=array([[-1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,-1.0]]); #inversion operator
    identity_mat=array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]); #identity operator
    identy_phase=exp(pi*2.j * dot(hkl, R)); #phase for initial site vector
    sym_phases=[];      #list of  symmetry operators (matrices) and the set of phases. Start with empty list.
    Rgen=rand(3);      #use random number to simulate general position to identify spacegroup forbidden reflections
    sum_phase_all=0;    #sum of all phases for site (to get scalar structure factor for site)
    sum_phase_gen=0;    #sum of phases for geneal (random) position
    for sym in spacegroup_list:
        mat=sym[0]
        vec=sym[1]
        time=sym[2]
        newR=dot(mat, R)+vec
        newRgen=dot(mat, Rgen)+vec
        phase=exp(pi*2.j * dot(hkl, newR))/identy_phase    #change phases so first one is unity
        #phaseexp(pi*2.j * dot(hkl, newR)); print 'phase: ',  phase;  print    #temp##################

        sum_phase_all+exp(pi*2.j * dot(hkl, newR))              #add new phase for site to sum
        sum_phase_gen+exp(pi*2.j * dot(hkl, newRgen))              #add new phase for general (random) position to sum
        newsym=1;
        for sym_phase in sym_phases:
            if allclose(mat,sym_phase[0], atol=tol) and allclose(time,sym_phase[1], atol=tol):              #compare mat with sym op is sym_phase list
                newsym=0;                                                       #if already in list then not new
                sym_phase[2]+=[phase];                                  #add new phase to phse list for sym op
                break
        if newsym==1:                                                         #sym op not in list so make a new entry
            sym_phases+=[[mat, time, [phase]]]


    sum_all_phases=0;                       #running total of all phases
    sum_Pplus_Tplus=0   #P even, T even etc
    sum_Pminus_Tplus=0
    sum_Pplus_Tminus=0
    sum_Pminus_Tminus=0
    Psym=Tsym=PTsym=None   #symmetries will be +1, -1 or 0 (even, odd, none)
    gen_scalar_allowed=site_scalar_allowed=1
    tensor_allowed=0
    if abs(sum_phase_all)<tol:
        site_scalar_allowed=0
    if abs(sum_phase_gen)<tol:
        gen_scalar_allowed=0        
    

    sum_phases=[]   #sum of phases for each symmetry operator
    for sym_phase in sym_phases:
        sum_phases+=[sum(sym_phase[2])]
        sum_all_phases+=sum(sym_phase[2])                             #add all phases (if all zero then forbidden for scalar)
        #if not allclose(sum(sym_phase[1]), 0, atol=tol):
        if not allclose(sum(sym_phase[2]), 0, atol=tol): ########## fix bug - hangover from before T was added
            tensor_allowed=1
        if allclose(sym_phase[0], identity_mat, atol=tol) and abs(sym_phase[1]-1)<tol:
            sum_Pplus_Tplus+=sum(sym_phase[2])
        elif allclose(sym_phase[0], inv_mat, atol=tol) and abs(sym_phase[1]-1)<tol:
            sum_Pminus_Tplus+=sum(sym_phase[2])
        elif allclose(sym_phase[0], identity_mat, atol=tol) and abs(sym_phase[1]+1)<tol: #time odd
            sum_Pplus_Tminus+=sum(sym_phase[2])                           
        elif allclose(sym_phase[0], inv_mat, atol=tol) and abs(sym_phase[1]+1)<tol:
            sum_Pminus_Tminus+=sum(sym_phase[2])                            

    print '<><><><><><><> sum_Pplus_Tminus <><><><>', sum_Pplus_Tminus

#    if tensor_allowed and abs(sum_Pplus_Tplus)>tol: #if there is no item with plus time and plus parity then there is no specific symmetry
    if tensor_allowed: #if there is no item with plus time and plus parity then there is no specific symmetry
        if sum_Pplus_Tplus-sum_Pminus_Tplus==0:
            Psym=+1
        if sum_Pplus_Tplus+sum_Pminus_Tplus==0:
            Psym=-1
        if sum_Pplus_Tplus-sum_Pplus_Tminus==0:
            Tsym=+1
        if sum_Pplus_Tplus+sum_Pplus_Tminus==0:
            Tsym=-1
        if sum_Pplus_Tplus-sum_Pminus_Tminus==0:
            PTsym=+1
        if sum_Pplus_Tplus+sum_Pminus_Tminus==0:
            PTsym=-1

    sum_phases=array(sum_phases)
    if allclose(sum_phases, sum_phases.real, atol=tol):
        sum_phases=real(sum_phases)
    else:
        print'=== Warning: sum of phases is compex. This was notexpected (see below):\n',sum_phases
    if abs(sum_phases[0])>tol:
        sum_phases=array(sum_phases); sum_phases=sum_phases/sum_phases[0] #normalize to first (identity)
    else:
        print '=== Warning: the phase sum for first (identity) operator is close to zero. This was notexpected'
            
    txtyn=['Yes','Invalid value', 'No', 'Invalid value']; txtoe=['Even', 'Odd', 'Either', 'Either']; 
    print '=== Analysis of structure actor symmetry phases\n','=== Sym_phases:',sym_phases, \
        '\n=== Site allowed:', msg(site_scalar_allowed, txt=txtyn), \
        '\n=== Spacegroup allowed:', msg(gen_scalar_allowed, txt=txtyn), \
        '\n=== Tensor allowed:', msg(tensor_allowed, txt=txtyn), \
        '\n=== Parity:', msg(Psym, txt=txtoe),'\n=== Time:', msg(Tsym, txt=txtoe),'\n=== PT:', msg(PTsym, txt=txtoe),\
        '\n=== sum_sym_phases:', sum_phases,\
        '\n=== Total ops:', len(sum_phases),\
        '\n=== Total positive:', len(sum_phases[sum_phases>0])

    sym_sum_phases=deepcopy(sym_phases)
    for ii in range(len(sym_sum_phases)):
        sym_sum_phases[ii][2]=sum_phases[ii]

    return [sym_sum_phases, sum_phases, gen_scalar_allowed, site_scalar_allowed, tensor_allowed, Psym, Tsym, PTsym, sym_phases,]

def msg(num, txt=['plus','minus','zero','other']):
    #return message text for +1,-1, 0, other (e.g. None)
    if num==1:
        str=txt[0]
    elif num==-1:
        str=txt[1]
    elif num==0:
        str=txt[2]
    else:
        str=txt[3]
    return str

def SF_symmetry_old(R, hkl, spacegroup_list):
    #analyse symmetry of any possible structure factor (mainly for information)
    inv_mat=array([[-1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,0.0,-1.0]]); #inversion operator
    identity_mat=array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]); #identity operator
    identy_phase=exp(pi*2.j * dot(hkl, R)); #phase for initial site vector
    sym_phases=[];      #list of  symmetry operators (matrices) and the set of phases. Start with empty list.
    Rgen=rand(3);      #use random number to simulate general position to identify spacegroup forbidden reflections
    sum_phase_all=0;    #sum of all phases for site (to get scalar structure factor for site)
    sum_phase_gen=0;    #sum of phases for geneal (random) position
    for sym in spacegroup_list:
        mat=sym[0]
        vec=sym[1]
        newR=dot(mat, R)+vec
        newRgen=dot(mat, Rgen)+vec
        phase=exp(pi*2.j * dot(hkl, newR))/identy_phase    #change phases so first one is unity
        #phase=exp(pi*2.j * dot(hkl, newR)); print 'phase: ',  phase;  print    #temp##################

        sum_phase_all+=exp(pi*2.j * dot(hkl, newR))              #add new phase for general (random) position to sum
        sum_phase_gen+=exp(pi*2.j * dot(hkl, newRgen))              #add new phase for general (random) position to sum
        newsym=1;
        for sym_phase in sym_phases:
            if allclose(mat,sym_phase[0], rtol=.001):              #compare mat with sym op is sym_phase list
                newsym=0;                                                       #if already in list then not new
                sym_phase[1]+=[phase];                                  #add new phase to phse list for sym op
                break
        if newsym==1:                                                         #sym op not in list so make a new entry
            sym_phases+=[[mat, [phase]]]
    #now analysis symmetries
    Podd=0;     #assumed not parity odd
    Peven=0;    #assumed not parity even
    Allowed=1;  #allowed untill proved otherwise
    sum_all_phases=0;                       #running total of all phases
    sum_identity_phases=0;               #running total of phases for indentity operator
    sum_inversion_phases=0;             #running total of phases for inversion operator
    all_phase_sums_zero=1;
    for sym_phase in sym_phases:
        sum_all_phases+=sum(sym_phase[1])                             #add all phases (if all zero then forbidden for scalar)
        if not allclose(sum(sym_phase[1]), 0, rtol=.001):
            all_phase_sums_zero=0;
        if allclose(sym_phase[0], identity_mat, rtol=.001):
            sum_identity_phases+=sum(sym_phase[1])                #add phases if  sym op is identity
        elif allclose(sym_phase[0], inv_mat, rtol=.001):
            sum_inversion_phases+=sum(sym_phase[1])               #add phases if  sym op is inversion

    scalar_forbidden=all_forbidden=parity_even=parity_odd=0;
    if allclose(sum_phase_gen, 0.0, rtol=.000001):
        sg_forbidden=1
        print '=== Spacegroup forbidden reflection (for scalar scattering) '
    else:
        sg_forbidden=0
        print '=== Not a spacegroup forbidden reflection (for scalar scattering) '
    if all_phase_sums_zero ==1:
        all_forbidden=1
        print '=== Forbidden for all tensors'
    else:
        all_forbidden=0
        print '=== Not forbidden for all tensors'
    if allclose(sum_phase_all, 0,  rtol=.0001):
        scalar_forbidden=1
        print '=== Forbidden reflection for specified sites (for scalar scattering)'
    else:
        print '=== Allowed reflection (for scalar scattering). Magnitude of site scalar structure factor:  %.5f' % abs(sum_phase_all)
        print sum_phase_all
    if allclose(sum_identity_phases-sum_inversion_phases, 0,  rtol=.001):
        parity_even=1
        print '=== Tensors are parity even'
    if allclose(sum_identity_phases+sum_inversion_phases, 0,  rtol=.001):
        parity_odd=1
        print '=== Tensors are parity odd'
    if parity_even==0 and parity_odd==0:
        print '=== Tensors may be parity odd or even'
    #print '=== Scalar structure factor for site: %6.2f' % sum_phase_all
    return [sym_phases, sg_forbidden,  scalar_forbidden,  all_forbidden,  parity_even,  parity_odd]

def latt2b(lat, direct=False, BLstyle=False):
    #follow Busing&Levy, D.E.Sands
    #direct=False: normal recip space B matrix (B&L)
    #direct=True, BLstyle=True: Busing & Levy style applied to real space (i.e. x||a)
    #direct=True, BLstyle=False: Real space B matrix compatible with recip space B matrix
    a1=lat[0];    a2=lat[1];    a3=lat[2];
    alpha1=lat[3]*pi/180;    alpha2=lat[4]*pi/180;    alpha3=lat[5]*pi/180;
    v=a1*a2*a3*sqrt(1-cos(alpha1)**2-cos(alpha2)**2-cos(alpha3)**2+2*cos(alpha1)*cos(alpha2)*cos(alpha3))
    b1=a2*a3*sin(alpha1)/v;    b2=a3*a1*sin(alpha2)/v;    b3=a1*a2*sin(alpha3)/v
    beta1=arccos((cos(alpha2)*cos(alpha3)-cos(alpha1))/sin(alpha2)/sin(alpha3))
    beta2=arccos((cos(alpha1)*cos(alpha3)-cos(alpha2))/sin(alpha3)/sin(alpha1))
    beta3=arccos((cos(alpha1)*cos(alpha2)-cos(alpha3))/sin(alpha1)/sin(alpha2))
    #reciprocal space
    B=array([    [b1, b2*cos(beta3), b3*cos(beta2)],
    [0, b2*sin(beta3), -b3*sin(beta2)*cos(alpha1)],
    [0, 0, 1/a3], ])
    #real space: Busing & Levy style applied to real space (i.e. x||a)
    BD=array([    [a1, a2*cos(alpha3), a3*cos(alpha2)],
    [0, a2*sin(alpha3), -a3*sin(alpha2)*cos(beta1)],
    [0, 0, 1/b3], ])
    # Real space  B matrix consistent with recip space B matrix (useful of calculations involve real and reiprocal space)
    Bdd=inv(B.transpose())

    if not direct:
        return B  
    else:
        if BLstyle:
            return BD
        else:
            return Bdd        

def spacegroup_list_from_genpos_list(genposlist):
    sglist=[];
    for genpos in genposlist:
        #sglist+=[genpos2matvec(genpos)]
        sglist+=[genpos2matvec(genpos)+[1]] #add +1 to indicate time symmetry
    return sglist

def site_sym(spacegroup_list, sitevec):
    symlist=[];
    tol=1e-6;   #coordinates treated as indintical if within tol
    sitevec=(sitevec+tol)%1-tol;      #map into range 0<=x<1 using tolerance tol
    for sg in spacegroup_list:
        #newpos=dot(sg[0], sitevec)+sg[1]
        newpos=dot(sg[0], sitevec)+sg[1]    #new coordinates after applying symmetry operator
        newpos=(newpos+tol)%1-tol;      #map into range 0<=x<1 using tolerance tol
        if allclose(newpos, sitevec, atol=.001):    #spacegroup operator presenves position so it is a point group operator
            symlist+=[[sg[0],sg[2]]]                                   #add matrix and time part of sg op to pg but...
            for sym in symlist[0:-1]:
                if allclose(sym[0], sg[0], atol=.001) and abs(sym[1]-sg[2])<tol:     #... remove it again if already in list
                    symlist=symlist[0: -1]
                    break
    return symlist

def site_sym_old(spacegroup_list, sitevec):
    symlist=[];
    tol=.001;   #coordinates treated as indintical if within tol
    sitevec=(sitevec+tol)%1-tol;      #map into range 0<=x<1 using tolerance tol
    for sg in spacegroup_list:
        #newpos=dot(sg[0], sitevec)+sg[1]
        newpos=dot(sg[0], sitevec)+sg[1]    #new coordinates after applying symmetry operator
        newpos=(newpos+tol)%1-tol;      #map into range 0<=x<1 using tolerance tol
        if allclose(newpos, sitevec, rtol=.001):    #spacegroup operator presenves position so it is a point group operator
            symlist+=[sg[0]];                                   #add matrix part of sg op to pg but...
            for sym in symlist[0:-1]:
                if allclose(sym, sg[0], rtol=.001):     #... remove it again if already in list
                    symlist=symlist[0: -1]
                    break
    return symlist


def equiv_sites(spacegroup_list, sitevec):
    poslist=[sitevec];
    tol=1e-6;   #coordinates treated as indintical if within tol
    for sg in spacegroup_list:
        newpos=dot(sg[0], sitevec)+sg[1]    #new coordinates after applying symmetry operator
        newpos=(newpos+tol)%1-tol;      #map into range 0<=x<1 using tolerance tol
        poslist+=[newpos]            #add new position to list...
        for pos in poslist[0:-1]:
            if allclose(pos, newpos, atol=tol):    #...if position already in list
                poslist=poslist[0:-1]                                   #...remove it
    return poslist

def crystal_point_sym(spacegroup_list):
    symlist=[]
    tol=1e-6;   #coordinates treated as indintical if within tol
    for sg in spacegroup_list:
        symlist+=[[sg[0],sg[2]]]                         #add matrix  and time part of sg op to pg but...
        for sym in symlist[0:-1]:
            if allclose(sym[0], sg[0], atol=tol) and abs(sym[1]-sg[2])<tol:     #... remove it again if already in list
                symlist=symlist[0: -1]
                break
    return symlist

def crystal_point_sym_old(spacegroup_list):
    symlist=[]
    for sg in spacegroup_list:
        symlist+=[sg[0]];                                   #add matrix part of sg op to pg but...
        for sym in symlist[0:-1]:
            if allclose(sym, sg[0], rtol=.001):     #... remove it again if already in list
                symlist=symlist[0: -1]
                break
    return symlist


def StoneSphericalToCartConversionCoefs(K,Calc=True,k=-1j):
    #Condon&Shortley phase convention (k=-i in Stone's paper)
    #from FortranForm (No - CForm?) First List->array, del other lists,spaces, extra bracket around first level
    #If Calc==False then use these expressions from Mathematica, else calculate them numerically
    if not Calc:
        if K==0:
            C=array(1.0)
        elif K==1:
            C=array(((Complex(0,1)/Sqrt(2),1/Sqrt(2),0),(0,0,Complex(0,1)), (Complex(0,-1)/Sqrt(2),1/Sqrt(2),0)))
        elif K==2:
            C=array((((-0.5,Complex(0,0.5),0),(Complex(0,0.5),0.5,0),(0,0,0)),((0,0,-0.5),(0,0,Complex(0,0.5)),(-0.5,Complex(0,0.5),0)),((1/Sqrt(6),0,0),(0,1/Sqrt(6),0),(0,0,-Sqrt(0.6666666666666666))),((0,0,0.5),(0,0,Complex(0,0.5)),(0.5,Complex(0,0.5),0)),((-0.5,Complex(0,-0.5),0),(Complex(0,-0.5),0.5,0),(0,0,0))))
        elif K==3:
            C=array(((((Complex(0,-0.5)/Sqrt(2),-1/(2.*Sqrt(2)),0),(-1/(2.*Sqrt(2)),Complex(0,0.5)/Sqrt(2),0),(0,0,0)),((-1/(2.*Sqrt(2)),Complex(0,0.5)/Sqrt(2),0),(Complex(0,0.5)/Sqrt(2),1/(2.*Sqrt(2)),0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0))),(((0,0,Complex(0,-0.5)/Sqrt(3)),(0,0,-1/(2.*Sqrt(3))),(Complex(0,-0.5)/Sqrt(3),-1/(2.*Sqrt(3)),0)),((0,0,-1/(2.*Sqrt(3))),(0,0,Complex(0,0.5)/Sqrt(3)),(-1/(2.*Sqrt(3)),Complex(0,0.5)/Sqrt(3),0)),((Complex(0,-0.5)/Sqrt(3),-1/(2.*Sqrt(3)),0),(-1/(2.*Sqrt(3)),Complex(0,0.5)/Sqrt(3),0),(0,0,0))),(((Complex(0,0.5)*Sqrt(0.3),1/(2.*Sqrt(30)),0),(1/(2.*Sqrt(30)),Complex(0,0.5)/Sqrt(30),0),(0,0,Complex(0,-1)*Sqrt(0.13333333333333333))),((1/(2.*Sqrt(30)),Complex(0,0.5)/Sqrt(30),0),(Complex(0,0.5)/Sqrt(30),Sqrt(0.3)/2.,0),(0,0,-Sqrt(0.13333333333333333))),((0,0,Complex(0,-1)*Sqrt(0.13333333333333333)),(0,0,-Sqrt(0.13333333333333333)),(Complex(0,-1)*Sqrt(0.13333333333333333),-Sqrt(0.13333333333333333),0))),(((0,0,Complex(0,1)/Sqrt(10)),(0,0,0),(Complex(0,1)/Sqrt(10),0,0)),((0,0,0),(0,0,Complex(0,1)/Sqrt(10)),(0,Complex(0,1)/Sqrt(10),0)),((Complex(0,1)/Sqrt(10),0,0),(0,Complex(0,1)/Sqrt(10),0),(0,0,Complex(0,-1)*Sqrt(0.4)))),(((Complex(0,-0.5)*Sqrt(0.3),1/(2.*Sqrt(30)),0),(1/(2.*Sqrt(30)),Complex(0,-0.5)/Sqrt(30),0),(0,0,Complex(0,1)*Sqrt(0.13333333333333333))),((1/(2.*Sqrt(30)),Complex(0,-0.5)/Sqrt(30),0),(Complex(0,-0.5)/Sqrt(30),Sqrt(0.3)/2.,0),(0,0,-Sqrt(0.13333333333333333))),((0,0,Complex(0,1)*Sqrt(0.13333333333333333)),(0,0,-Sqrt(0.13333333333333333)),(Complex(0,1)*Sqrt(0.13333333333333333),-Sqrt(0.13333333333333333),0))),(((0,0,Complex(0,-0.5)/Sqrt(3)),(0,0,1/(2.*Sqrt(3))),(Complex(0,-0.5)/Sqrt(3),1/(2.*Sqrt(3)),0)),((0,0,1/(2.*Sqrt(3))),(0,0,Complex(0,0.5)/Sqrt(3)),(1/(2.*Sqrt(3)),Complex(0,0.5)/Sqrt(3),0)),((Complex(0,-0.5)/Sqrt(3),1/(2.*Sqrt(3)),0),(1/(2.*Sqrt(3)),Complex(0,0.5)/Sqrt(3),0),(0,0,0))),(((Complex(0,0.5)/Sqrt(2),-1/(2.*Sqrt(2)),0),(-1/(2.*Sqrt(2)),Complex(0,-0.5)/Sqrt(2),0),(0,0,0)),((-1/(2.*Sqrt(2)),Complex(0,-0.5)/Sqrt(2),0),(Complex(0,-0.5)/Sqrt(2),1/(2.*Sqrt(2)),0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0)))))
        elif K==4:
            C=array((((((0.25,Complex(0,-0.25),0),(Complex(0,-0.25),-0.25,0),(0,0,0)),((Complex(0,-0.25),-0.25,0),(-0.25,Complex(0,0.25),0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0))),(((Complex(0,-0.25),-0.25,0),(-0.25,Complex(0,0.25),0),(0,0,0)),((-0.25,Complex(0,0.25),0),(Complex(0,0.25),0.25,0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0))),(((0,0,0),(0,0,0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0)))),((((0,0,1/(4.*Sqrt(2))),(0,0,Complex(0,-0.25)/Sqrt(2)),(1/(4.*Sqrt(2)),Complex(0,-0.25)/Sqrt(2),0)),((0,0,Complex(0,-0.25)/Sqrt(2)),(0,0,-1/(4.*Sqrt(2))),(Complex(0,-0.25)/Sqrt(2),-1/(4.*Sqrt(2)),0)),((1/(4.*Sqrt(2)),Complex(0,-0.25)/Sqrt(2),0),(Complex(0,-0.25)/Sqrt(2),-1/(4.*Sqrt(2)),0),(0,0,0))),(((0,0,Complex(0,-0.25)/Sqrt(2)),(0,0,-1/(4.*Sqrt(2))),(Complex(0,-0.25)/Sqrt(2),-1/(4.*Sqrt(2)),0)),((0,0,-1/(4.*Sqrt(2))),(0,0,Complex(0,0.25)/Sqrt(2)),(-1/(4.*Sqrt(2)),Complex(0,0.25)/Sqrt(2),0)),((Complex(0,-0.25)/Sqrt(2),-1/(4.*Sqrt(2)),0),(-1/(4.*Sqrt(2)),Complex(0,0.25)/Sqrt(2),0),(0,0,0))),(((1/(4.*Sqrt(2)),Complex(0,-0.25)/Sqrt(2),0),(Complex(0,-0.25)/Sqrt(2),-1/(4.*Sqrt(2)),0),(0,0,0)),((Complex(0,-0.25)/Sqrt(2),-1/(4.*Sqrt(2)),0),(-1/(4.*Sqrt(2)),Complex(0,0.25)/Sqrt(2),0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0)))),((((-1/(2.*Sqrt(7)),Complex(0,0.25)/Sqrt(7),0),(Complex(0,0.25)/Sqrt(7),0,0),(0,0,1/(2.*Sqrt(7)))),((Complex(0,0.25)/Sqrt(7),0,0),(0,Complex(0,0.25)/Sqrt(7),0),(0,0,Complex(0,-0.5)/Sqrt(7))),((0,0,1/(2.*Sqrt(7))),(0,0,Complex(0,-0.5)/Sqrt(7)),(1/(2.*Sqrt(7)),Complex(0,-0.5)/Sqrt(7),0))),(((Complex(0,0.25)/Sqrt(7),0,0),(0,Complex(0,0.25)/Sqrt(7),0),(0,0,Complex(0,-0.5)/Sqrt(7))),((0,Complex(0,0.25)/Sqrt(7),0),(Complex(0,0.25)/Sqrt(7),1/(2.*Sqrt(7)),0),(0,0,-1/(2.*Sqrt(7)))),((0,0,Complex(0,-0.5)/Sqrt(7)),(0,0,-1/(2.*Sqrt(7))),(Complex(0,-0.5)/Sqrt(7),-1/(2.*Sqrt(7)),0))),(((0,0,1/(2.*Sqrt(7))),(0,0,Complex(0,-0.5)/Sqrt(7)),(1/(2.*Sqrt(7)),Complex(0,-0.5)/Sqrt(7),0)),((0,0,Complex(0,-0.5)/Sqrt(7)),(0,0,-1/(2.*Sqrt(7))),(Complex(0,-0.5)/Sqrt(7),-1/(2.*Sqrt(7)),0)),((1/(2.*Sqrt(7)),Complex(0,-0.5)/Sqrt(7),0),(Complex(0,-0.5)/Sqrt(7),-1/(2.*Sqrt(7)),0),(0,0,0)))),((((0,0,-3/(4.*Sqrt(14))),(0,0,Complex(0,0.25)/Sqrt(14)),(-3/(4.*Sqrt(14)),Complex(0,0.25)/Sqrt(14),0)),((0,0,Complex(0,0.25)/Sqrt(14)),(0,0,-1/(4.*Sqrt(14))),(Complex(0,0.25)/Sqrt(14),-1/(4.*Sqrt(14)),0)),((-3/(4.*Sqrt(14)),Complex(0,0.25)/Sqrt(14),0),(Complex(0,0.25)/Sqrt(14),-1/(4.*Sqrt(14)),0),(0,0,1/Sqrt(14)))),(((0,0,Complex(0,0.25)/Sqrt(14)),(0,0,-1/(4.*Sqrt(14))),(Complex(0,0.25)/Sqrt(14),-1/(4.*Sqrt(14)),0)),((0,0,-1/(4.*Sqrt(14))),(0,0,Complex(0,0.75)/Sqrt(14)),(-1/(4.*Sqrt(14)),Complex(0,0.75)/Sqrt(14),0)),((Complex(0,0.25)/Sqrt(14),-1/(4.*Sqrt(14)),0),(-1/(4.*Sqrt(14)),Complex(0,0.75)/Sqrt(14),0),(0,0,Complex(0,-1)/Sqrt(14)))),(((-3/(4.*Sqrt(14)),Complex(0,0.25)/Sqrt(14),0),(Complex(0,0.25)/Sqrt(14),-1/(4.*Sqrt(14)),0),(0,0,1/Sqrt(14))),((Complex(0,0.25)/Sqrt(14),-1/(4.*Sqrt(14)),0),(-1/(4.*Sqrt(14)),Complex(0,0.75)/Sqrt(14),0),(0,0,Complex(0,-1)/Sqrt(14))),((0,0,1/Sqrt(14)),(0,0,Complex(0,-1)/Sqrt(14)),(1/Sqrt(14),Complex(0,-1)/Sqrt(14),0)))),((((3/(2.*Sqrt(70)),0,0),(0,1/(2.*Sqrt(70)),0),(0,0,-Sqrt(0.05714285714285714))),((0,1/(2.*Sqrt(70)),0),(1/(2.*Sqrt(70)),0,0),(0,0,0)),((0,0,-Sqrt(0.05714285714285714)),(0,0,0),(-Sqrt(0.05714285714285714),0,0))),(((0,1/(2.*Sqrt(70)),0),(1/(2.*Sqrt(70)),0,0),(0,0,0)),((1/(2.*Sqrt(70)),0,0),(0,3/(2.*Sqrt(70)),0),(0,0,-Sqrt(0.05714285714285714))),((0,0,0),(0,0,-Sqrt(0.05714285714285714)),(0,-Sqrt(0.05714285714285714),0))),(((0,0,-Sqrt(0.05714285714285714)),(0,0,0),(-Sqrt(0.05714285714285714),0,0)),((0,0,0),(0,0,-Sqrt(0.05714285714285714)),(0,-Sqrt(0.05714285714285714),0)),((-Sqrt(0.05714285714285714),0,0),(0,-Sqrt(0.05714285714285714),0),(0,0,2*Sqrt(0.05714285714285714))))),((((0,0,3/(4.*Sqrt(14))),(0,0,Complex(0,0.25)/Sqrt(14)),(3/(4.*Sqrt(14)),Complex(0,0.25)/Sqrt(14),0)),((0,0,Complex(0,0.25)/Sqrt(14)),(0,0,1/(4.*Sqrt(14))),(Complex(0,0.25)/Sqrt(14),1/(4.*Sqrt(14)),0)),((3/(4.*Sqrt(14)),Complex(0,0.25)/Sqrt(14),0),(Complex(0,0.25)/Sqrt(14),1/(4.*Sqrt(14)),0),(0,0,-(1/Sqrt(14))))),(((0,0,Complex(0,0.25)/Sqrt(14)),(0,0,1/(4.*Sqrt(14))),(Complex(0,0.25)/Sqrt(14),1/(4.*Sqrt(14)),0)),((0,0,1/(4.*Sqrt(14))),(0,0,Complex(0,0.75)/Sqrt(14)),(1/(4.*Sqrt(14)),Complex(0,0.75)/Sqrt(14),0)),((Complex(0,0.25)/Sqrt(14),1/(4.*Sqrt(14)),0),(1/(4.*Sqrt(14)),Complex(0,0.75)/Sqrt(14),0),(0,0,Complex(0,-1)/Sqrt(14)))),(((3/(4.*Sqrt(14)),Complex(0,0.25)/Sqrt(14),0),(Complex(0,0.25)/Sqrt(14),1/(4.*Sqrt(14)),0),(0,0,-(1/Sqrt(14)))),((Complex(0,0.25)/Sqrt(14),1/(4.*Sqrt(14)),0),(1/(4.*Sqrt(14)),Complex(0,0.75)/Sqrt(14),0),(0,0,Complex(0,-1)/Sqrt(14))),((0,0,-(1/Sqrt(14))),(0,0,Complex(0,-1)/Sqrt(14)),(-(1/Sqrt(14)),Complex(0,-1)/Sqrt(14),0)))),((((-1/(2.*Sqrt(7)),Complex(0,-0.25)/Sqrt(7),0),(Complex(0,-0.25)/Sqrt(7),0,0),(0,0,1/(2.*Sqrt(7)))),((Complex(0,-0.25)/Sqrt(7),0,0),(0,Complex(0,-0.25)/Sqrt(7),0),(0,0,Complex(0,0.5)/Sqrt(7))),((0,0,1/(2.*Sqrt(7))),(0,0,Complex(0,0.5)/Sqrt(7)),(1/(2.*Sqrt(7)),Complex(0,0.5)/Sqrt(7),0))),(((Complex(0,-0.25)/Sqrt(7),0,0),(0,Complex(0,-0.25)/Sqrt(7),0),(0,0,Complex(0,0.5)/Sqrt(7))),((0,Complex(0,-0.25)/Sqrt(7),0),(Complex(0,-0.25)/Sqrt(7),1/(2.*Sqrt(7)),0),(0,0,-1/(2.*Sqrt(7)))),((0,0,Complex(0,0.5)/Sqrt(7)),(0,0,-1/(2.*Sqrt(7))),(Complex(0,0.5)/Sqrt(7),-1/(2.*Sqrt(7)),0))),(((0,0,1/(2.*Sqrt(7))),(0,0,Complex(0,0.5)/Sqrt(7)),(1/(2.*Sqrt(7)),Complex(0,0.5)/Sqrt(7),0)),((0,0,Complex(0,0.5)/Sqrt(7)),(0,0,-1/(2.*Sqrt(7))),(Complex(0,0.5)/Sqrt(7),-1/(2.*Sqrt(7)),0)),((1/(2.*Sqrt(7)),Complex(0,0.5)/Sqrt(7),0),(Complex(0,0.5)/Sqrt(7),-1/(2.*Sqrt(7)),0),(0,0,0)))),((((0,0,-1/(4.*Sqrt(2))),(0,0,Complex(0,-0.25)/Sqrt(2)),(-1/(4.*Sqrt(2)),Complex(0,-0.25)/Sqrt(2),0)),((0,0,Complex(0,-0.25)/Sqrt(2)),(0,0,1/(4.*Sqrt(2))),(Complex(0,-0.25)/Sqrt(2),1/(4.*Sqrt(2)),0)),((-1/(4.*Sqrt(2)),Complex(0,-0.25)/Sqrt(2),0),(Complex(0,-0.25)/Sqrt(2),1/(4.*Sqrt(2)),0),(0,0,0))),(((0,0,Complex(0,-0.25)/Sqrt(2)),(0,0,1/(4.*Sqrt(2))),(Complex(0,-0.25)/Sqrt(2),1/(4.*Sqrt(2)),0)),((0,0,1/(4.*Sqrt(2))),(0,0,Complex(0,0.25)/Sqrt(2)),(1/(4.*Sqrt(2)),Complex(0,0.25)/Sqrt(2),0)),((Complex(0,-0.25)/Sqrt(2),1/(4.*Sqrt(2)),0),(1/(4.*Sqrt(2)),Complex(0,0.25)/Sqrt(2),0),(0,0,0))),(((-1/(4.*Sqrt(2)),Complex(0,-0.25)/Sqrt(2),0),(Complex(0,-0.25)/Sqrt(2),1/(4.*Sqrt(2)),0),(0,0,0)),((Complex(0,-0.25)/Sqrt(2),1/(4.*Sqrt(2)),0),(1/(4.*Sqrt(2)),Complex(0,0.25)/Sqrt(2),0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0)))),((((0.25,Complex(0,0.25),0),(Complex(0,0.25),-0.25,0),(0,0,0)),((Complex(0,0.25),-0.25,0),(-0.25,Complex(0,-0.25),0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0))),(((Complex(0,0.25),-0.25,0),(-0.25,Complex(0,-0.25),0),(0,0,0)),((-0.25,Complex(0,-0.25),0),(Complex(0,-0.25),0.25,0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0))),(((0,0,0),(0,0,0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0)),((0,0,0),(0,0,0),(0,0,0))))))
        else:
            raise ValueError('No Spherical to Cart conversion availble for rank '+str(K))
    else:
        CS=[];
        for i in range(1,K+1):      #generate coupling sequence CS=[1,2,3...K]
            CS+=[i]
        C=array(StoneCoefficients(CS,k=k)).transpose()   
    return C

def spherical_to_cart_tensor(Ts):
    K=(len(Ts)-1)/2; #spherical tensor rank
    C=StoneSphericalToCartConversionCoefs(K).conjugate()
    Tc=C[0]*0.0;    #array of zeros
    for kk in range(-K, K+1):
        Tc=Tc+Ts[kk+K]*C[kk+K]
    return Tc

def cart_to_spherical_tensor(Tc):
    K=len(Tc.shape); #Cartesian tensor rank
    Cconj=StoneSphericalToCartConversionCoefs(K)
    Ts=zeros(2*K+1, dtype=complex)
    for kk in range(-K, K+1):
        Ts[kk+K]=sum(Cconj[kk+K]*Tc)
        #print kk+K, Ts, Ts[kk+K],  sum(Cconj[kk+K]*Tc)
    return Ts

def theta_to_cartesian(hkl,hkln,psi,B):
    Rx=array([
    [1, 0, 0],
    [0, cos(psi), sin(psi)],
    [0, -sin(psi), cos(psi)]
    ])
    #cross only seems to work with row vectors - take transpose

    xp=dot(B, hkl)/norm(dot(B, hkl))
    cp=cross(dot(B, hkl),dot(B, hkln));
    zp=cp/norm(cp)
    yp=cross(zp, xp);
    Ucpsi=array([xp, yp, zp]).T;
    #Ucpsi=array([xp, yp, zp])##testing -  not right
    Uctheta=dot(Ucpsi, Rx)#original
    return Uctheta

def Complex(a, b):          #allow FortranForm Mathematica output with only minor mods need to delete spaces, List( ... )->array([...])
    return a+b*1j

def Sqrt(a):
    return sqrt(a)

def scalar_contract(X, T):
    if len(X)!=len(T):
        raise ValueError("Can't form scalar contraction of tensors with different rank")
    K=len(X)/2
    scalar=0
    for kk in range(len(X)):
        q=kk-K
        scalar+=(-1)**q * X[2*K-kk] * T[kk]
    return scalar

def Xtensor(process, rank, time, parity, e0, e1, q0, q1):
    [e0x, e0y, e0z]=e0
    [e1x, e1y, e1z]=e1
    [q0x, q0y, q0z]=q0
    [q1x, q1y, q1z]=q1

    if process=='E1E1' and rank==2 and parity==1:
        X2=array([
        ((e0x-Complex(0,1)*e0y)*(e1x-Complex(0,1)*e1y))/2.,
        (e0z*(e1x-Complex(0,1)*e1y)+(e0x-Complex(0,1)*e0y)*e1z)/2.,
        -((e0x*e1x+e0y*e1y-2*e0z*e1z)/Sqrt(6)),
        (-(e0z*e1x)-Complex(0,1)*e0z*e1y-e0x*e1z-Complex(0,1)*e0y*e1z)/2.,
        ((e0x+Complex(0,1)*e0y)*(e1x+Complex(0,1)*e1y))/2.
        ])
        return X2
    if process=='E1E2' and rank==3 and parity==-1:
        #E1E2 rank 3
        n3t=array([
        ((e0x-Complex(0,1)*e0y)*(Complex(0,1)*e1x+e1y)*
        (q0x-Complex(0,1)*q0y))/(2.*Sqrt(10)),
        (e0z*(Complex(0,1)*e1x+e1y)*(q0x-Complex(0,1)*q0y)+
        (Complex(0,1)*e0x+e0y)*
        (e1z*(q0x-Complex(0,1)*q0y)+(e1x-Complex(0,1)*e1y)*q0z))/
        (2.*Sqrt(15)),(Complex(0,0.1)*
        (-((e0x-Complex(0,1)*e0y)*(e1x+Complex(0,1)*e1y)*
        (q0x-Complex(0,1)*q0y))+
        4*e1z*(e0z*(q0x-Complex(0,1)*q0y)+
        (e0x-Complex(0,1)*e0y)*q0z)-
        2*(e1x-Complex(0,1)*e1y)*(e0x*q0x+e0y*q0y-2*e0z*q0z)))/Sqrt(6),
        (Complex(0,-0.2)*(e0x*e1z*q0x+e0y*e1z*q0y+e0x*e1x*q0z+
        e0y*e1y*q0z+e0z*(e1x*q0x+e1y*q0y-2*e1z*q0z)))/Sqrt(2),
        (Complex(0,0.1)*((e0x+Complex(0,1)*e0y)*(e1x-Complex(0,1)*e1y)*
        (q0x+Complex(0,1)*q0y)-
        4*e1z*(e0z*(q0x+Complex(0,1)*q0y)+
        (e0x+Complex(0,1)*e0y)*q0z)+
        2*(e1x+Complex(0,1)*e1y)*(e0x*q0x+e0y*q0y-2*e0z*q0z)))/Sqrt(6),
        (Complex(0,0.5)*(e0z*(e1x+Complex(0,1)*e1y)*(q0x+Complex(0,1)*q0y)+
        (e0x+Complex(0,1)*e0y)*
        (e1z*(q0x+Complex(0,1)*q0y)+(e1x+Complex(0,1)*e1y)*q0z)))/
        Sqrt(15),((e0x+Complex(0,1)*e0y)*(Complex(0,-1)*e1x+e1y)*
        (q0x+Complex(0,1)*q0y))/(2.*Sqrt(10))
        ])
        n3=array([
        ((e0x-Complex(0,1)*e0y)*(Complex(0,1)*e1x+e1y)*
        (q1x-Complex(0,1)*q1y))/(2.*Sqrt(10)),
        (e0z*(Complex(0,1)*e1x+e1y)*(q1x-Complex(0,1)*q1y)+
        (Complex(0,1)*e0x+e0y)*
        (e1z*(q1x-Complex(0,1)*q1y)+(e1x-Complex(0,1)*e1y)*q1z))/
        (2.*Sqrt(15)),(Complex(0,0.1)*
        (-((e0x+Complex(0,1)*e0y)*(e1x-Complex(0,1)*e1y)*
        (q1x-Complex(0,1)*q1y))+
        4*e0z*(e1z*(q1x-Complex(0,1)*q1y)+
        (e1x-Complex(0,1)*e1y)*q1z)-
        2*(e0x-Complex(0,1)*e0y)*(e1x*q1x+e1y*q1y-2*e1z*q1z)))/Sqrt(6),
        (Complex(0,-0.2)*(e0x*e1z*q1x+e0y*e1z*q1y+e0x*e1x*q1z+
        e0y*e1y*q1z+e0z*(e1x*q1x+e1y*q1y-2*e1z*q1z)))/Sqrt(2),
        (Complex(0,0.1)*((e0x-Complex(0,1)*e0y)*(e1x+Complex(0,1)*e1y)*
        (q1x+Complex(0,1)*q1y)-
        4*e0z*(e1z*(q1x+Complex(0,1)*q1y)+
        (e1x+Complex(0,1)*e1y)*q1z)+
        2*(e0x+Complex(0,1)*e0y)*(e1x*q1x+e1y*q1y-2*e1z*q1z)))/Sqrt(6),
        (Complex(0,0.5)*(e0z*(e1x+Complex(0,1)*e1y)*(q1x+Complex(0,1)*q1y)+
        (e0x+Complex(0,1)*e0y)*
        (e1z*(q1x+Complex(0,1)*q1y)+(e1x+Complex(0,1)*e1y)*q1z)))/
        Sqrt(15),((e0x+Complex(0,1)*e0y)*(Complex(0,-1)*e1x+e1y)*
        (q1x+Complex(0,1)*q1y))/(2.*Sqrt(10))
        ])
        if time==1:
            return n3t-n3
        elif time==-1:
            return  n3t+n3
    if process=='E1E2' and rank==2 and parity==-1:
        #E1E2 rank 2
        n2t=array([(e0z*(Complex(0,1)*e1x+e1y)*(q0x-Complex(0,1)*q0y)+(e0x-Complex(0,1)*e0y)*(Complex(0,-2)*e1z*q0x-2*e1z*q0y+Complex(0,1)*e1x*q0z+e1y*q0z))/(2.*Sqrt(30)),(e0z*(Complex(0,-1)*e1z*q0x-e1z*q0y+Complex(0,2)*e1x*q0z+2*e1y*q0z)+e0y*(Complex(0,1)*e1y*q0x+e1x*(q0x-Complex(0,2)*q0y)-e1z*q0z)+e0x*(-2*e1y*q0x+e1x*q0y+Complex(0,1)*e1y*q0y-Complex(0,1)*e1z*q0z))/(2.*Sqrt(30)),(-(e0z*e1y*q0x)+e0z*e1x*q0y+e0y*e1x*q0z-e0x*e1y*q0z)/(2.*Sqrt(5)),(e0z*(Complex(0,-1)*e1z*q0x+e1z*q0y+Complex(0,2)*e1x*q0z-2*e1y*q0z)+e0x*(2*e1y*q0x-e1x*q0y+Complex(0,1)*e1y*q0y-Complex(0,1)*e1z*q0z)+e0y*(Complex(0,1)*e1y*q0x-e1x*(q0x+Complex(0,2)*q0y)+e1z*q0z))/(2.*Sqrt(30)),(e0z*(Complex(0,-1)*e1x+e1y)*(q0x+Complex(0,1)*q0y)+(e0x+Complex(0,1)*e0y)*(Complex(0,2)*e1z*q0x-2*e1z*q0y-Complex(0,1)*e1x*q0z+e1y*q0z))/(2.*Sqrt(30))])
        n2=array([(Complex(0,-2)*e0z*(e1x-Complex(0,1)*e1y)*(q1x-Complex(0,1)*q1y)+(Complex(0,1)*e0x+e0y)*(e1z*(q1x-Complex(0,1)*q1y)+(e1x-Complex(0,1)*e1y)*q1z))/(2.*Sqrt(30)),(Complex(0,-1)*e0z*(e1z*q1x-Complex(0,1)*e1z*q1y+e1x*q1z-Complex(0,1)*e1y*q1z)+e0x*(e1y*q1x+e1x*q1y-Complex(0,2)*e1y*q1y+Complex(0,2)*e1z*q1z)+e0y*(-2*e1x*q1x+Complex(0,1)*e1y*q1x+Complex(0,1)*e1x*q1y+2*e1z*q1z))/(2.*Sqrt(30)),(-(e0y*(e1z*q1x+e1x*q1z))+e0x*(e1z*q1y+e1y*q1z))/(2.*Sqrt(5)),(e0z*(Complex(0,-1)*e1z*q1x+e1z*q1y-Complex(0,1)*e1x*q1z+e1y*q1z)+e0y*(2*e1x*q1x+Complex(0,1)*e1y*q1x+Complex(0,1)*e1x*q1y-2*e1z*q1z)-e0x*(e1y*q1x+e1x*q1y+Complex(0,2)*e1y*q1y-Complex(0,2)*e1z*q1z))/(2.*Sqrt(30)),(Complex(0,2)*e0z*(e1x+Complex(0,1)*e1y)*(q1x+Complex(0,1)*q1y)+(Complex(0,-1)*e0x+e0y)*(e1z*(q1x+Complex(0,1)*q1y)+(e1x+Complex(0,1)*e1y)*q1z))/(2.*Sqrt(30))])
        if time==1:
            return n2t-n2
        elif time==-1:
            return  n2t+n2
    if process=='E1E2' and rank==1 and parity==-1:
        #E1E2 rank 1
        print "xxxxx E1E2 rank 1 not tested"
        n1t=array([(Complex(0,0.1)*(-3*(e0x-Complex(0,1)*e0y)*(e1x+Complex(0,1)*e1y)*(q0x-Complex(0,1)*q0y)-3*e1z*(e0z*(q0x-Complex(0,1)*q0y)+(e0x-Complex(0,1)*e0y)*q0z)-(e1x-Complex(0,1)*e1y)*(e0x*q0x+e0y*q0y-2*e0z*q0z)))/Sqrt(6),(Complex(0,-0.1)*(-2*e0x*e1z*q0x-2*e0y*e1z*q0y+3*e0x*e1x*q0z+3*e0y*e1y*q0z+e0z*(3*e1x*q0x+3*e1y*q0y+4*e1z*q0z)))/Sqrt(3),(Complex(0,0.1)*(3*(e0x+Complex(0,1)*e0y)*(e1x-Complex(0,1)*e1y)*(q0x+Complex(0,1)*q0y)+3*e1z*(e0z*(q0x+Complex(0,1)*q0y)+(e0x+Complex(0,1)*e0y)*q0z)+(e1x+Complex(0,1)*e1y)*(e0x*q0x+e0y*q0y-2*e0z*q0z)))/Sqrt(6)])
        n1=array([(Complex(0,0.1)*(-3*(e0x+Complex(0,1)*e0y)*(e1x-Complex(0,1)*e1y)*(q1x-Complex(0,1)*q1y)-3*e0z*(e1z*(q1x-Complex(0,1)*q1y)+(e1x-Complex(0,1)*e1y)*q1z)-(e0x-Complex(0,1)*e0y)*(e1x*q1x+e1y*q1y-2*e1z*q1z)))/Sqrt(6),(Complex(0,0.1)*(-3*(e0x*e1z*q1x+e0y*e1z*q1y+e0x*e1x*q1z+e0y*e1y*q1z)+2*e0z*(e1x*q1x+e1y*q1y-2*e1z*q1z)))/Sqrt(3),(Complex(0,0.1)*(3*(e0x-Complex(0,1)*e0y)*(e1x+Complex(0,1)*e1y)*(q1x+Complex(0,1)*q1y)+3*e0z*(e1z*(q1x+Complex(0,1)*q1y)+(e1x+Complex(0,1)*e1y)*q1z)+(e0x+Complex(0,1)*e0y)*(e1x*q1x+e1y*q1y-2*e1z*q1z)))/Sqrt(6)])
        if time==1:
            return n1t-n1
        elif time==-1:
            return  n1t+n1

    if process=='E2E2' and rank==0 and parity==1 and time==1:
        XQQ0=array([(e0y*(3*e1y*q0x*q1x-2*e1x*q0y*q1x+3*e1x*q0x*q1y+4*e1y*q0y*q1y+
3*e1z*q0z*q1y-2*e1z*q0y*q1z+3*e1y*q0z*q1z)+
e0z*(3*e1z*q0x*q1x-2*e1x*q0z*q1x+3*e1z*q0y*q1y-2*e1y*q0z*q1y+
3*e1x*q0x*q1z+3*e1y*q0y*q1z+4*e1z*q0z*q1z)+
e0x*(3*e1y*q0y*q1x+3*e1z*q0z*q1x-2*e1y*q0x*q1y-2*e1z*q0x*q1z+
e1x*(4*q0x*q1x+3*q0y*q1y+3*q0z*q1z)))/(6.*Sqrt(5))])
        return XQQ0
    if process=='E2E2' and rank==1 and parity==1 and time==-1:
        XQQ1=array([(e0x*(-2*e1x*q0z*q1x+Complex(0,1)*e1y*q0z*q1x+
Complex(0,1)*e1x*q0z*q1y+2*e1x*q0x*q1z-
Complex(0,1)*e1x*q0y*q1z+e1y*q0y*q1z+
e1z*(2*q0x*q1x-Complex(0,1)*q0y*q1x+q0y*q1y+2*q0z*q1z))-
Complex(0,1)*e0y*(Complex(0,-1)*e1y*q0z*q1x-
Complex(0,1)*e1x*q0z*q1y-2*e1y*q0z*q1y+e1x*q0x*q1z+
Complex(0,1)*e1y*q0x*q1z+2*e1y*q0y*q1z+
e1z*(q0x*q1x+Complex(0,1)*q0x*q1y+2*q0y*q1y+2*q0z*q1z))-
e0z*(e1x*(2*q0x*q1x-Complex(0,1)*q0x*q1y+q0y*q1y+2*q0z*q1z)-
Complex(0,1)*(2*e1z*(Complex(0,1)*q0z*q1x+q0z*q1y-
Complex(0,1)*q0x*q1z-q0y*q1z)+
e1y*(q0x*q1x+Complex(0,1)*q0y*q1x+2*q0y*q1y+2*q0z*q1z))))/
(4.*Sqrt(5)),(Complex(0,0.5)*
(e0z*(e1z*q0y*q1x-e1z*q0x*q1y-e1y*q0x*q1z+e1x*q0y*q1z)-
e0x*(2*e1y*q0x*q1x-2*e1x*q0y*q1x+2*e1x*q0x*q1y+
2*e1y*q0y*q1y+e1z*q0z*q1y+e1y*q0z*q1z)+
e0y*(2*e1y*q0y*q1x+e1z*q0z*q1x-2*e1y*q0x*q1y+
e1x*(2*q0x*q1x+2*q0y*q1y+q0z*q1z))))/Sqrt(10),
(e0x*(e1x*(-2*q0z*q1x-Complex(0,1)*q0z*q1y+2*q0x*q1z+
Complex(0,1)*q0y*q1z)+e1y*(Complex(0,-1)*q0z*q1x+q0y*q1z)+
e1z*(2*q0x*q1x+Complex(0,1)*q0y*q1x+q0y*q1y+2*q0z*q1z))-
e0z*(2*e1z*(q0z*(q1x+Complex(0,1)*q1y)-
(q0x+Complex(0,1)*q0y)*q1z)+
e1y*(Complex(0,1)*q0x*q1x+q0y*q1x+Complex(0,2)*q0y*q1y+
Complex(0,2)*q0z*q1z)+
e1x*(2*q0x*q1x+Complex(0,1)*q0x*q1y+q0y*q1y+2*q0z*q1z))+
e0y*(-(e1x*q0z*q1y)+Complex(0,1)*e1x*q0x*q1z+
e1y*(-(q0z*(q1x+Complex(0,2)*q1y))+
(q0x+Complex(0,2)*q0y)*q1z)+
e1z*(q0x*(Complex(0,1)*q1x+q1y)+
Complex(0,2)*(q0y*q1y+q0z*q1z))))/(4.*Sqrt(5))])
        return XQQ1
    if process=='E2E2' and rank==2 and parity==1 and time==1:
        XQQ2=array([(e0z*(-3*e1z*(q0x-Complex(0,1)*q0y)*(q1x-Complex(0,1)*q1y)+
(e1x-Complex(0,1)*e1y)*
(4*q0z*q1x-Complex(0,4)*q0z*q1y-3*q0x*q1z+
Complex(0,3)*q0y*q1z))+
e0y*(e1z*(Complex(0,3)*q0z*q1x+3*q0z*q1y-Complex(0,4)*q0x*q1z-
4*q0y*q1z)+Complex(0,1)*e1x*
(2*q0x*q1x+2*q0y*q1y+3*q0z*q1z)+
e1y*(Complex(0,2)*q0y*q1x+Complex(0,2)*q0x*q1y+4*q0y*q1y+
3*q0z*q1z))+e0x*(e1x*
(-4*q0x*q1x+Complex(0,2)*q0y*q1x+Complex(0,2)*q0x*q1y-
3*q0z*q1z)+Complex(0,1)*
(e1z*(Complex(0,3)*q0z*q1x+3*q0z*q1y-Complex(0,4)*q0x*q1z-
4*q0y*q1z)+e1y*(2*q0x*q1x+2*q0y*q1y+3*q0z*q1z))))/
(4.*Sqrt(21)),(e0y*(e1y*(-3*q0z*q1x+Complex(0,2)*q0z*q1y-
3*q0x*q1z+Complex(0,2)*q0y*q1z)+
e1x*(Complex(0,-4)*q0z*q1x-3*q0z*q1y+Complex(0,3)*q0x*q1z+
4*q0y*q1z)+e1z*(Complex(0,3)*q0x*q1x+4*q0y*q1x-
3*q0x*q1y+Complex(0,2)*q0y*q1y+Complex(0,2)*q0z*q1z))-
e0x*(e1x*(2*q0z*q1x-Complex(0,3)*q0z*q1y+2*q0x*q1z-
Complex(0,3)*q0y*q1z)+
e1y*(Complex(0,-3)*q0z*q1x-4*q0z*q1y+Complex(0,4)*q0x*q1z+
3*q0y*q1z)+e1z*(Complex(0,-3)*q0y*q1x+
2*q0x*(q1x+Complex(0,2)*q1y)+3*q0y*q1y+2*q0z*q1z))+
e0z*(-2*e1z*(q0z*q1x-Complex(0,1)*q0z*q1y+q0x*q1z-
Complex(0,1)*q0y*q1z)+
e1y*(Complex(0,3)*q0x*q1x-3*q0y*q1x+4*q0x*q1y+
Complex(0,2)*q0y*q1y+Complex(0,2)*q0z*q1z)-
e1x*(2*q0x*q1x+Complex(0,4)*q0y*q1x-Complex(0,3)*q0x*q1y+
3*q0y*q1y+2*q0z*q1z)))/(4.*Sqrt(21)),
(e0y*(6*e1y*q0x*q1x-8*e1x*q0y*q1x+6*e1x*q0x*q1y+4*e1y*q0y*q1y-
3*e1z*q0z*q1y+4*e1z*q0y*q1z-3*e1y*q0z*q1z)-
e0z*(3*e1z*q0x*q1x-4*e1x*q0z*q1x+3*e1z*q0y*q1y-4*e1y*q0z*q1y+
3*e1x*q0x*q1z+3*e1y*q0y*q1z+8*e1z*q0z*q1z)+
e0x*(6*e1y*q0y*q1x-3*e1z*q0z*q1x-8*e1y*q0x*q1y+4*e1z*q0x*q1z+
e1x*(4*q0x*q1x+6*q0y*q1y-3*q0z*q1z)))/(6.*Sqrt(14)),
(e0y*(e1x*(Complex(0,-4)*q0z*q1x+3*q0z*q1y+Complex(0,3)*q0x*q1z-
4*q0y*q1z)+e1y*(3*q0z*q1x+Complex(0,2)*q0z*q1y+
3*q0x*q1z+Complex(0,2)*q0y*q1z)+
e1z*(-4*q0y*q1x+Complex(0,2)*q0y*q1y+
3*q0x*(Complex(0,1)*q1x+q1y)+Complex(0,2)*q0z*q1z))+
e0x*(e1x*(2*q0z*q1x+Complex(0,3)*q0z*q1y+2*q0x*q1z+
Complex(0,3)*q0y*q1z)+
e1y*(Complex(0,3)*q0z*q1x-4*q0z*q1y-Complex(0,4)*q0x*q1z+
3*q0y*q1z)+e1z*(Complex(0,3)*q0y*q1x+
2*q0x*(q1x-Complex(0,2)*q1y)+3*q0y*q1y+2*q0z*q1z))+
e0z*(2*e1z*(q0z*q1x+Complex(0,1)*q0z*q1y+q0x*q1z+
Complex(0,1)*q0y*q1z)+
e1y*(Complex(0,3)*q0x*q1x+3*q0y*q1x-4*q0x*q1y+
Complex(0,2)*q0y*q1y+Complex(0,2)*q0z*q1z)+
e1x*(2*q0x*q1x-Complex(0,4)*q0y*q1x+Complex(0,3)*q0x*q1y+
3*q0y*q1y+2*q0z*q1z)))/(4.*Sqrt(21)),
(Complex(0,-1)*(e0z*(3*e1z*(Complex(0,-1)*q0x+q0y)*
(q1x+Complex(0,1)*q1y)+
(e1x+Complex(0,1)*e1y)*
(Complex(0,4)*q0z*q1x-4*q0z*q1y-Complex(0,3)*q0x*q1z+
3*q0y*q1z))+e0y*
(e1z*(3*q0z*(q1x+Complex(0,1)*q1y)-
4*(q0x+Complex(0,1)*q0y)*q1z)+
e1y*(2*q0y*(q1x+Complex(0,2)*q1y)+2*q0x*q1y+
Complex(0,3)*q0z*q1z)+
e1x*(2*q0x*q1x+2*q0y*q1y+3*q0z*q1z)))-
e0x*(e1x*(4*q0x*q1x+Complex(0,2)*q0y*q1x+Complex(0,2)*q0x*q1y+
3*q0z*q1z)+Complex(0,1)*
(e1z*(Complex(0,-3)*q0z*q1x+3*q0z*q1y+Complex(0,4)*q0x*q1z-
4*q0y*q1z)+e1y*(2*q0x*q1x+2*q0y*q1y+3*q0z*q1z))))/
(4.*Sqrt(21))])
        return XQQ2
    if process=='E2E2' and rank==3 and parity==1 and time==-1:
        XQQ3=array([(e0z*(e1x-Complex(0,1)*e1y)*(q0x-Complex(0,1)*q0y)*
(q1x-Complex(0,1)*q1y)-
(e0x-Complex(0,1)*e0y)*e1z*(q0x-Complex(0,1)*q0y)*
(q1x-Complex(0,1)*q1y)+
(e0x-Complex(0,1)*e0y)*(e1x-Complex(0,1)*e1y)*q0z*
(q1x-Complex(0,1)*q1y)-
(e0x-Complex(0,1)*e0y)*(e1x-Complex(0,1)*e1y)*
(q0x-Complex(0,1)*q0y)*q1z)/(4.*Sqrt(2)),
(Complex(0,0.25)*(Complex(0,-2)*e0z*(e1x-Complex(0,1)*e1y)*q0z*
(q1x-Complex(0,1)*q1y)+
e0y*(e1y*q0y*q1x-e1y*q0x*q1y+
e1x*(-(q0x*q1x)+Complex(0,2)*q0y*q1x+q0y*q1y)+
2*e1z*q0x*q1z-Complex(0,2)*e1z*q0y*q1z)+
e0x*(e1y*q0x*q1x-e1x*q0y*q1x+e1x*q0x*q1y-
Complex(0,2)*e1y*q0x*q1y-e1y*q0y*q1y+
Complex(0,2)*e1z*q0x*q1z+2*e1z*q0y*q1z)))/Sqrt(3),
(Complex(0,-1)*e0y*(e1x*(5*q0z*q1x-Complex(0,3)*q0z*q1y+
3*q0x*q1z-Complex(0,5)*q0y*q1z)+
e1y*(Complex(0,-3)*q0z*q1x-q0z*q1y+Complex(0,3)*q0x*q1z+
q0y*q1z)+e1z*(Complex(0,-5)*q0y*q1x+
3*q0x*(q1x+Complex(0,1)*q1y)+q0y*q1y-4*q0z*q1z))+
e0x*(e1x*(-(q0z*q1x)+Complex(0,3)*q0z*q1y+q0x*q1z-
Complex(0,3)*q0y*q1z)+
e1y*(Complex(0,3)*q0z*q1x+5*q0z*q1y+Complex(0,5)*q0x*q1z+
3*q0y*q1z)+e1z*(Complex(0,-3)*q0y*q1x+
q0x*(q1x+Complex(0,5)*q1y)+3*q0y*q1y-4*q0z*q1z))+
e0z*(4*e1z*(q0z*q1x-Complex(0,1)*q0z*q1y-q0x*q1z+
Complex(0,1)*q0y*q1z)-
e1x*(Complex(0,5)*q0y*q1x+q0x*(q1x-Complex(0,3)*q1y)+
3*q0y*q1y-4*q0z*q1z)+
e1y*(Complex(0,3)*q0x*q1x-3*q0y*q1x+5*q0x*q1y+
Complex(0,1)*q0y*q1y-Complex(0,4)*q0z*q1z)))/(4.*Sqrt(30)),
(Complex(0,0.5)*(2*e0z*(-(e1z*q0y*q1x)+e1z*q0x*q1y+e1y*q0x*q1z-
e1x*q0y*q1z)+e0y*(e1y*q0y*q1x-2*e1z*q0z*q1x-e1y*q0x*q1y+
e1x*(q0x*q1x+q0y*q1y-2*q0z*q1z))+
e0x*(e1x*q0y*q1x-e1x*q0x*q1y+2*e1z*q0z*q1y-
e1y*(q0x*q1x+q0y*q1y-2*q0z*q1z))))/Sqrt(10),
(Complex(0,1)*e0y*(e1x*(5*q0z*q1x+Complex(0,3)*q0z*q1y+3*q0x*q1z+
Complex(0,5)*q0y*q1z)+
e1y*(Complex(0,3)*q0z*q1x-q0z*q1y-Complex(0,3)*q0x*q1z+
q0y*q1z)+e1z*(Complex(0,5)*q0y*q1x+
3*q0x*(q1x-Complex(0,1)*q1y)+q0y*q1y-4*q0z*q1z))+
e0x*(e1x*(-(q0z*(q1x+Complex(0,3)*q1y))+
(q0x+Complex(0,3)*q0y)*q1z)+
e1y*(Complex(0,-3)*q0z*q1x+5*q0z*q1y-Complex(0,5)*q0x*q1z+
3*q0y*q1z)+e1z*(Complex(0,3)*q0y*q1x+
q0x*(q1x-Complex(0,5)*q1y)+3*q0y*q1y-4*q0z*q1z))-
e0z*(4*e1z*(-(q0z*(q1x+Complex(0,1)*q1y))+
(q0x+Complex(0,1)*q0y)*q1z)+
e1x*(Complex(0,-5)*q0y*q1x+q0x*(q1x+Complex(0,3)*q1y)+
3*q0y*q1y-4*q0z*q1z)+
e1y*(Complex(0,3)*q0x*q1x+3*q0y*q1x-5*q0x*q1y+
Complex(0,1)*q0y*q1y-Complex(0,4)*q0z*q1z)))/(4.*Sqrt(30)),
(Complex(0,0.25)*(Complex(0,2)*e0z*(e1x+Complex(0,1)*e1y)*q0z*
(q1x+Complex(0,1)*q1y)+
e0y*(e1y*q0y*q1x-e1y*q0x*q1y+
e1x*(-(q0x*q1x)-Complex(0,2)*q0y*q1x+q0y*q1y)+
2*e1z*q0x*q1z+Complex(0,2)*e1z*q0y*q1z)+
e0x*(e1y*q0x*q1x-e1x*q0y*q1x+e1x*q0x*q1y+
Complex(0,2)*e1y*q0x*q1y-e1y*q0y*q1y-
Complex(0,2)*e1z*q0x*q1z+2*e1z*q0y*q1z)))/Sqrt(3),
(e0z*(e1x+Complex(0,1)*e1y)*(q0x+Complex(0,1)*q0y)*
(q1x+Complex(0,1)*q1y)-
(e0x+Complex(0,1)*e0y)*(e1z*(q0x+Complex(0,1)*q0y)*
(q1x+Complex(0,1)*q1y)-
(e1x+Complex(0,1)*e1y)*
(q0z*(q1x+Complex(0,1)*q1y)-(q0x+Complex(0,1)*q0y)*q1z)))/
(4.*Sqrt(2))])
        return XQQ3
    if process=='E2E2' and rank==4 and parity==1 and time==1:
        #nedit replace \n \r and ' ' with nothing.
        #List((( -> array([(( and add square bracket after first round, also at end
        XQQ4=array([((e0x-Complex(0,1)*e0y)*(e1x-Complex(0,1)*e1y)*
(q0x-Complex(0,1)*q0y)*(q1x-Complex(0,1)*q1y))/4.,
(e0z*(e1x-Complex(0,1)*e1y)*(q0x-Complex(0,1)*q0y)*
(q1x-Complex(0,1)*q1y)+
(e0x-Complex(0,1)*e0y)*e1z*(q0x-Complex(0,1)*q0y)*
(q1x-Complex(0,1)*q1y)+
(e0x-Complex(0,1)*e0y)*(e1x-Complex(0,1)*e1y)*q0z*
(q1x-Complex(0,1)*q1y)+
(e0x-Complex(0,1)*e0y)*(e1x-Complex(0,1)*e1y)*
(q0x-Complex(0,1)*q0y)*q1z)/(4.*Sqrt(2)),
(e0x*(e1x*(-2*q0x*q1x+Complex(0,1)*q0y*q1x+Complex(0,1)*q0x*q1y+
2*q0z*q1z)+Complex(0,1)*
(Complex(0,-2)*e1z*(q0z*q1x-Complex(0,1)*q0z*q1y+q0x*q1z-
Complex(0,1)*q0y*q1z)+e1y*(q0x*q1x+q0y*q1y-2*q0z*q1z))\
)+Complex(0,1)*(Complex(0,-2)*e0z*
(e1z*(q0x-Complex(0,1)*q0y)*(q1x-Complex(0,1)*q1y)+
(e1x-Complex(0,1)*e1y)*
(q0z*(q1x-Complex(0,1)*q1y)+(q0x-Complex(0,1)*q0y)*q1z))\
+e0y*(-2*e1z*(q0z*q1x-Complex(0,1)*q0z*q1y+q0x*q1z-
Complex(0,1)*q0y*q1z)+
e1x*(q0x*q1x+q0y*q1y-2*q0z*q1z)+
e1y*(q0y*q1x+q0x*q1y-Complex(0,2)*q0y*q1y+
Complex(0,2)*q0z*q1z))))/(4.*Sqrt(7)),
(-(e0x*(e1x*(3*q0z*q1x-Complex(0,1)*q0z*q1y+3*q0x*q1z-
Complex(0,1)*q0y*q1z)+
e1y*(Complex(0,-1)*q0z*q1x+q0z*q1y-Complex(0,1)*q0x*q1z+
q0y*q1z)+e1z*(3*q0x*q1x-Complex(0,1)*q0y*q1x-
Complex(0,1)*q0x*q1y+q0y*q1y-4*q0z*q1z)))-
e0z*(-4*e1z*(q0z*q1x-Complex(0,1)*q0z*q1y+q0x*q1z-
Complex(0,1)*q0y*q1z)+
e1x*(3*q0x*q1x-Complex(0,1)*q0y*q1x-Complex(0,1)*q0x*q1y+
q0y*q1y-4*q0z*q1z)+
e1y*(Complex(0,-1)*q0x*q1x+q0y*q1x+q0x*q1y-
Complex(0,3)*q0y*q1y+Complex(0,4)*q0z*q1z))-
e0y*(e1y*(q0z*q1x-Complex(0,3)*q0z*q1y+q0x*q1z-
Complex(0,3)*q0y*q1z)+
e1x*(Complex(0,-1)*q0z*q1x+q0z*q1y-Complex(0,1)*q0x*q1z+
q0y*q1z)+e1z*(q0y*q1x-Complex(0,3)*q0y*q1y+
q0x*(Complex(0,-1)*q1x+q1y)+Complex(0,4)*q0z*q1z)))/
(4.*Sqrt(14)),(e0y*(e1y*q0x*q1x+e1x*q0y*q1x+e1x*q0x*q1y+
3*e1y*q0y*q1y-4*e1z*q0z*q1y-4*e1z*q0y*q1z-4*e1y*q0z*q1z)-
4*e0z*(e1z*q0x*q1x+e1x*q0z*q1x+e1z*q0y*q1y+e1y*q0z*q1y+
e1x*q0x*q1z+e1y*q0y*q1z-2*e1z*q0z*q1z)+
e0x*(e1y*q0y*q1x-4*e1z*q0z*q1x+e1y*q0x*q1y-4*e1z*q0x*q1z+
e1x*(3*q0x*q1x+q0y*q1y-4*q0z*q1z)))/(2.*Sqrt(70)),
(e0x*(e1x*(3*q0z*q1x+Complex(0,1)*q0z*q1y+3*q0x*q1z+
Complex(0,1)*q0y*q1z)+
e1y*(Complex(0,1)*q0z*q1x+q0z*q1y+Complex(0,1)*q0x*q1z+
q0y*q1z)+e1z*(3*q0x*q1x+Complex(0,1)*q0y*q1x+
Complex(0,1)*q0x*q1y+q0y*q1y-4*q0z*q1z))+
e0z*(-4*e1z*(q0z*q1x+Complex(0,1)*q0z*q1y+q0x*q1z+
Complex(0,1)*q0y*q1z)+
e1x*(3*q0x*q1x+Complex(0,1)*q0y*q1x+Complex(0,1)*q0x*q1y+
q0y*q1y-4*q0z*q1z)+
e1y*(Complex(0,1)*q0x*q1x+q0y*q1x+q0x*q1y+
Complex(0,3)*q0y*q1y-Complex(0,4)*q0z*q1z))+
e0y*(e1y*(q0z*q1x+Complex(0,3)*q0z*q1y+q0x*q1z+
Complex(0,3)*q0y*q1z)+
e1x*(Complex(0,1)*q0z*q1x+q0z*q1y+Complex(0,1)*q0x*q1z+
q0y*q1z)+e1z*(q0y*q1x+Complex(0,3)*q0y*q1y+
q0x*(Complex(0,1)*q1x+q1y)-Complex(0,4)*q0z*q1z)))/
(4.*Sqrt(14)),(2*e0z*(e1z*(q0x+Complex(0,1)*q0y)*
(q1x+Complex(0,1)*q1y)+
(e1x+Complex(0,1)*e1y)*
(q0z*(q1x+Complex(0,1)*q1y)+(q0x+Complex(0,1)*q0y)*q1z))-
Complex(0,1)*e0y*(-2*e1z*
(q0z*q1x+Complex(0,1)*q0z*q1y+q0x*q1z+
Complex(0,1)*q0y*q1z)+e1x*(q0x*q1x+q0y*q1y-2*q0z*q1z)+
e1y*(q0y*q1x+q0x*q1y+Complex(0,2)*q0y*q1y-
Complex(0,2)*q0z*q1z))+
e0x*(-(e1x*(2*q0x*q1x+Complex(0,1)*q0y*q1x+
Complex(0,1)*q0x*q1y-2*q0z*q1z))-
Complex(0,1)*(Complex(0,2)*e1z*
(q0z*q1x+Complex(0,1)*q0z*q1y+q0x*q1z+
Complex(0,1)*q0y*q1z)+e1y*(q0x*q1x+q0y*q1y-2*q0z*q1z)))\
)/(4.*Sqrt(7)),(-(e0z*(e1x+Complex(0,1)*e1y)*(q0x+Complex(0,1)*q0y)*
(q1x+Complex(0,1)*q1y))-
(e0x+Complex(0,1)*e0y)*
(e1z*(q0x+Complex(0,1)*q0y)*(q1x+Complex(0,1)*q1y)+
(e1x+Complex(0,1)*e1y)*
(q0z*(q1x+Complex(0,1)*q1y)+(q0x+Complex(0,1)*q0y)*q1z)))/
(4.*Sqrt(2)),((e0x+Complex(0,1)*e0y)*(e1x+Complex(0,1)*e1y)*
(q0x+Complex(0,1)*q0y)*(q1x+Complex(0,1)*q1y))/4.])
        return XQQ4
    else:
        raise ValueError('Unknown tensor type')

def norm_array(Array, Minval=0.001):
    #Normalise array by largest abs value if >Minval (avoids trying to renormalise zero array)
    greatest=0.0
    flatarray=Array.flat
    for i in range(len(flatarray)):
        if abs(flatarray[i])>abs(greatest):
            greatest=flatarray[i]
    if abs(greatest)>Minval:
        newarray=Array/greatest
    else:
        newarray=Array
    return newarray

def array_interp(arraylist,n):
    '''
    Interpolates arraylist with n extra vectors linearly spaced between each entry; converts lists to arrays
    '''
    outlist=[]
    for ii in range(len(arraylist)-1):
        inc=(array(arraylist[ii+1])-array(arraylist[ii]))/(n+1.0)
        for jj in range(n+1):
            outlist+=[array(arraylist[ii])+inc*jj]
    outlist+=[array(arraylist[-1])]
    return outlist

def magvector_fourier(fourier_coeffs,site_vector):
    '''
    calculate magnetic vector at particular site based on a list of vector fourier components
    modulation vectors are in reciprocal lattice units (crystal coordinate system)
    magnetic vectors are in Cartesian coordinates
    each member of the 'fourier_coeffs' list is a two element list comprising:
        - a real moldulation vector (array)
        - a complex amplitide vector (array)
    '''
    magvec=zeros(3,dtype=complex)
    for modulation in fourier_coeffs:
        tau=modulation[0]
        amp=modulation[1]
        #print tau, amp, exp(2*pi*1.j*dot(tau,site_vector))*amp
        magvec+=exp(2*pi*1.j*dot(tau,site_vector))*amp
    return magvec


def calc_magSF_fourier(fourier_coeffs,sitevec,spacegroup,hkl):
    '''
    calculate magnetic unit cell structure factor for hkl using magvector_fourier function
    uses spacegroup and sitevec to calculate all sites in unit cell (single equivalent ion)
    '''
    tol=1e-6
    SF=0
    sites=[]
    for sym in spacegroup:
        mat=sym[0]
        vec=sym[1]
        newR=dot(mat, sitevec)+vec
        newR=(newR+tol)%1-tol;      #map into range 0<=x<1 using tolerance tol
        newsite=True
        for site in sites:
            if allclose(newR, site, atol=tol):    #site is already in list
                newsite=False
                break
        if newsite:    
            phase=exp(pi*2.j * dot(hkl, newR))
            magvec=magvector_fourier(fourier_coeffs,newR)
            print "site, dot(hkl, newR), phase, magvec: ", newR, dot(hkl, newR), phase, magvec
            SF+=phase*magvec
            sites+=[newR]
            #print sites
    return SF
    
## routines adapted from Sypy code    
    
def ClebschGordan(j1, j2, m1, m2, J, M, warn=True):
    """
    ClebschGordan(j1, j2, m1, m2, J, M, cglimit=20,warn=True)
    Computes exact sympy form for Clebsch-Gordan coefficient
    <j1 j2; m1 m2|j1 j2; JM>.
    For reference see
    http://en.wikipedia.org/wiki/Table_of_Clebsch-Gordan_coefficients.
    Clebsch Gordan numpy function by Michael V. DePalatis, modified and converted to sympy by SPC
    warn gives warning for unphysical coefficients
     Adapted from sympy ClebschGordan
    """
    j1=float(j1); j2=float(j2); m1=float(m1); m2=float(m2); J=float(J); M=float(M);
    if not M==(m1+m2) or J>(j1+j2) or J<abs(j1-j2) or J<0 or abs(m1)>j1 or abs(m2)>j2 or abs(M)>J:
        if warn:
            print 'Warning: Unphysical Clebsch-Gordan coefficient (j1,j2,m1,m2,J,M)='+str((j1,j2,m1,m2,J,M))
        return 0

    c1 = sqrt((2*J+1) * factorial(J+j1-j2) * factorial(J-j1+j2) * \
        factorial(j1+j2-J)/factorial(j1+j2+J+1))
    c2 = sqrt(factorial(J+M) * factorial(J-M) * factorial(j1-m1) * \
        factorial(j1+m1) * factorial(j2-m2) * factorial(j2+m2))
    c3 = 0.
    cglimit=max((j1+j2-J),(j1-m1),(j2+m1))+1        #max k that satisfies requirement that all factorial args are non-neg
    for k in arange(cglimit):
        use = True
        d = [0, 0, 0, 0, 0]
        d[0] = j1 + j2 - J - k
        d[1] = j1 - m1 - k
        d[2] = j2 + m2 - k
        d[3] = J - j2 + m1 + k
        d[4] = J - j1 -m2 + k
        prod = factorial(k)
        for arg in d:
            if arg < 0:
                use = False
                break
            prod *= factorial(arg)
        if use:
            #print k
            c3 += (-1)**k/prod
    return c1*c2*c3


def StoneCoefficients(CouplingSequenceList,k=-1j):
    '''
    StoneCoefficients(CouplingSequenceList,k=phase_convention)
    Sympy Spherical-Cartesian conversion coefficients from
    A.J. Stone Molecular Physics 29 1461 (1975) (Equation 1.9)
    CouplingSequenceList is the coupling sequence for spherical tensors, 
        each time coupling to a new vector to form a tensor of given rank
        (sequence always starts with 1)
    k=-I for Condon & Shortley phase convention (default) of k=1 for Racah
    e.g. StoneCoefficients([1,2,3]) returns conversion coefficients for K=3, coupling with 
    maximum rank and Condon & Shortley (default) phase convention
    Example:     C123=StoneCoefficients([1,2,3])    returns conversion matrix for coupling sequence 123 (K=3)
            print lcontract(C123,3,[1,0,0,0,0,0,0]) returns table values for Q=-3
    Numpy version converted from, Sympy version
    '''
    rt=2**0.5               #sqrt(2)
    #Cartesian index sequence: x,y,z; spherical index sequence -q...+q
    C1=[[1j*k/rt,0,-1j*k/rt],[k/rt,0,k/rt],[0,1j*k,0]]    #coefficients for vector (K=1)=C1
    N=len(CouplingSequenceList)                #total number of vectors coupled    
    #diag('line 311',['rt','C1','N'],locals()) 
    if N==0:
        C=[1]    
    elif N>0:
        C=C1
    if N>1:
        if CouplingSequenceList[0]!=1:
            raise ValueError('First rank in sequence must be 1')
        for J in CouplingSequenceList[1:]:        #loop through all J's after the first
            Cnew=StoneCoupleVector(C,J,C1)        #couple to next vector to make tensor of rank J
            C=Cnew     
            #diag('stone coef main loop',['C','J','C1','Cnew'],locals()) 
    return C


def StoneCoupleVector(Cold,Knew,C1):
    '''
    StoneCoupleVector(Cold,Knew,C1)
    couple Stone coefficients Cold to a new vector to make coefficient for spherocal tensor of rank Knew
    using vector coupling coefficients C1
    A.J. Stone Molecular Physics 29 1461 (1975) (Equation 1.9)
    Numpy version converted from Sympy version
    '''
    #indexing of Cartesian components=0,1,2 (x,y,z); indexing of spherical components=m+j (0..2j)
    Cold=array(Cold); C1=array(C1);                  #make sure they are arrays
    oldshape=Cold.shape                            #shape of previous connversion matrix
    newshape=len(oldshape)*[3]+[2*Knew+1]                    #shape of new conversion matrix
    Cnew=zeros(newshape, dtype=complex)                            #empty matrix for new conversion matrix
    oldindlist=indexlist(oldshape)                       #list of all indices for old matrix
    jn=(newshape[-1]-1)/2
    jn_=(oldshape[-1]-1)/2                            #j_{n-1} from Stone
    for ind in oldindlist:                            #loop through indices
        ind=list(ind);                  #make sure its a list
        mp=ind[-1]-jn_
        oldelement=Cold[tuple(ind)]                #value of old matrix for index
        for an in range(3):

            for mpp in [-1,0,1]:
                m=mp+mpp                    #definitions follow Stone...
                vectorconversion=C1[an][mpp+1]            #required element of vector conversion matrix 
                cj=ClebschGordan(jn_,1,mp,mpp,jn,m)
                element=Cold[tuple(ind)]        #required element of old conversion matrix 
                newind=ind[:-1]+[an]+[m+jn]            #index for new matrix
                newelement=Cnew[tuple(newind)]    #element of new matrix
                newelement+=element*vectorconversion*cj        #add in new bit
                Cnew[tuple(newind)]=copy(newelement)    #save back to array
                #diag_tmp= element*vectorconversion*cj               
                #diag('stonecouplevector main loop',['newind','cj','element','vectorconversion','diag_tmp','newelement','m','cj'],locals()) 


    return Cnew

    
def indexlist(shape):
	'''
	indexlist(shape)
	create a list of index lists covering all indices for shape list (all possible indices)
      (Numpy 1.6 has new indexing functionality that my render this obsolete)
	'''
	lenshape=len(shape)				#number of indices
	totalelements=1					#total number of elements in nested list (start value)
	indlistlist=[]					#start with empty list of index lists
	indlist=list(range(lenshape))			#single index list placeholder
	multiplicity=list(range(lenshape))		#multiplicities of collumns for counting 
	multiplicity[lenshape-1]=1			#last collumn is 'ones'	
	for i in range(lenshape-2,-1,-1):		
		totalelements*=shape[i+1]
		multiplicity[i]=multiplicity[i+1]*shape[i+1]
	totalelements*=shape[0]

	for count in range(totalelements):		#loop through each element (flat index)
		for i in range(lenshape-1,-1,-1):	#loop through indices in reverse order
			indlist[i]=(count/multiplicity[i])%shape[i]	
		indlistlist+=[copy(indlist)]
	return indlistlist    
    
#def getNestedListItem(listarray,listindex):
#    '''
#    getNestedListItem(listarray,listindex)
#    Return element of nested list listarray given by indices in listindex
#    '''
#    indexstring=''
#    for i in listindex:
#        indexstring+='['+str(int(i))+']'
#    return eval('listarray'+indexstring)
#
#def setNestedListItem(listarray,listindex,element):
#    '''
#    setNestedListItem(listarray,listindex,element)
#    Set element of nested list listarray given by indices in listindex to element
#    '''
#    indexstring=''
#    for i in listindex:
#        indexstring+='['+str(int(i))+']'
#    exec('listarray'+indexstring+'=element')
#
#def lshape(listarray):
#	'''
#	lshape(listarray)
#	returns a list giving the shape of a multi-dimensional nested list (multidimensional array)
#	only the first elements are looked at - there is no check that the shape is consistent
#	returns empty list if listarray has no len method (e.g. a scalar)
#	'''
#	shape=[]
#	while True:
#		try:
#			shape+=[len(listarray)]
#			listarray=listarray[0]
#		except:
#			return shape
#   
#def lzeros(shapelist,number=0):
#	'''
#	lzeros(shapelist,number=0)
#	creates nested lists of zeros (or specified number, symbol etc) of shape given by shapelist
#	'''
#
#	newlist=list(shapelist)
#	newlist.reverse()			#reverse a copy of the shape list
#	m=number
#	for dim in newlist:
#		newm=[]
#		for i in range(dim):
#			newm+=[copy(m)]
#		m=newm
#	return m
# 

#
##for test
#StoneCoefficients([1,2])
##set break point at offending line
##c (continiue) about 14 times till error
