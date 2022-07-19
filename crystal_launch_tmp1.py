import Crystal
import os

class xtal(Crystal.Crystal):
    def mercury(self):
        print 'Launch external application mercury with crytal data'
        tmp_cif_file='/tmp/crystal_tmp_cif_file.cif'
        self.save_cif(tmp_cif_file)
        os.system('module load mercury; mercury '+tmp_cif_file)


cy=xtal()
#cy.load_cif('/home/spc93/spc_cifs/CoCO3_icsd_61066.cif')
cy.load_cif('/home/spc93/spc_cifs/BiFeO3_r3c_icsd_168319.cif')
#cy.load_cif('/home/spc93/spc_cifs/Fe3BO6 icsd_1910.cif')
#cy.load_cif('/home/spc93/spc_cifs/GGG_sawada_icsd_84874.cif')


cy.reflection_list(5)
cy.mercury()

# most cif files fail

