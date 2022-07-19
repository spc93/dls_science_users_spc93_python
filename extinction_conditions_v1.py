'''
Does the lattice matter for extinction conditions or should it be chosen to be consistent with spacegroup?
'''



import TensorScatteringClass as ten
import numpy as np


class ExtinctionConditions(ten.TensorScatteringClass):  
    '''
    Extend TensorScatteringClass to look at extinction conditions
    Start with normal spacegroup condidtions      
    '''

    pass  

def isRat(r):
    return r in [0, 1, 1/2, 1/3, 2/3, 1/4, 3/4]

t = ExtinctionConditions(spacegroup_number = 18, wyckoff_letter = 'a', lattice = [1.0, 1.0, 1.0, 90, 90, 90])
sglist = t.spacegroup_list_from_genpos_list(t.symxyz) 
print(sglist)

R_wyck = t.sitevec
R_gen = np.array([np.random.rand(), np.random.rand(), np.random.rand()])
hkl = np.array([2, 1, 0])
print('\nhkl: ', hkl)
print('Wyckoff position', R_wyck)
print(t.SF_symmetry(R_wyck, hkl, sglist))
print('\nGeneral position', R_gen)
print(t.SF_symmetry(R_gen, hkl, sglist))

print(t.allR)

n_atoms = len(t.allR)



smallestNonRat = 1
for i in range(n_atoms):
    for j in range(n_atoms):
        z = t.allR[j][0] - t.allR[i][0]
        if not isRat(z-np.floor(z)):
            if z-np.floor(z) < smallestNonRat:
                smallestNonRat = z-np.floor(z)
        print(i, j,  z, z-np.floor(z), isRat(z-np.floor(z)), smallestNonRat)

   

R, r = np.zeros(n_atoms), np.zeros(n_atoms)
 
print('smallestNonRat:', smallestNonRat)
   
    

def findrat(x1, x2):
    #if x2 can be written R +/- x1 then write R and sign
    if isRat(x2-x1):
        print('rational part:', x2-x1, 'sign:', +1)
    elif isRat(x2+x1):
        print('rational part:', x2+x1, 'sign:', -1)
    else:
        print('No relationship')



findrat(smallestNonRat, 1/2-smallestNonRat)



############## either create symbolic coords from floats or do the while thing symbolically...........


#[i/j for i in range(3) for j in range(1, 3)]

#function to test for simple rational

# then find smallest non-rational

#then write all coords as ratinal +/- nonrational

# write each atomic coord as rational + a, where a is irrational