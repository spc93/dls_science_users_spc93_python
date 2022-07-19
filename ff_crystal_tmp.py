import Crystal
import MagneticSpacegroup

CoCO3=MagneticSpacegroup.MagneticSpacegroup()

CoCO3.load_cif('/dls_sw/i16/software/python/userscripts/steve/subgroup_cif.txt',axes_from_parent_cell=np.transpose(np.array([[-2,-1,0],[0,1,0],[2/3,1/3,1/3]])), origin_from_parent_cell=[-1/6,-1/3,1/6])

#Calculate the reflections:
#CoCO3.reflection_list(5.223)
#Cfr. Charge / parent / magnetic reflections
#Subtraction and addiction of reflection lists