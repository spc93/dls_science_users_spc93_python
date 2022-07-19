
# write equivalent symInfoFromCIF
# use both in tensor calcs

def symInfoFromSpacegroupAndWyckoff(sg_num, wyck_letter, site = None):
    '''
    Get symxyz (general positions as strings from spacegroup and Wyckoff letter
    returns site after substituting either given x, y, z values or random numbers
    Requires CCTBX package
    '''
    import cctbx.sgtbx as sg

    sgi = sg.space_group_info(sg_num)
    sg = sgi.group()
    w = sgi.wyckoff_table()
    
    all_letters = list([w.position(i).letter() for i in range(w.size())])
    
    if wyck_letter in all_letters:
        opxyz = str(w.position(wyck_letter).special_op()) # generic xyz string for site
    else:
        print('=== Invalid Wyckoff letter for spacegroup. Valid letters:', all_letters)
        raise(ValueError)

    opxyz = opxyz.replace('/','./') # fix for Python 2.x int devide - not needed for Python 3

    if not site == None:        #use random values for x,y,z if not specified
        x, y, z = site
    else:
        x, y, z = rand(), rand(), rand()

    sitevec = np.array(eval(opxyz)) # numerical array for site
    symxyz = [str(s.as_xyz()) for s in sg.all_ops()]    # all sg ops as xyz list
    sitestr = '%s %s' % (str(w.position(wyck_letter).multiplicity()), wyck_letter) 

    return (symxyz, sitevec, sitestr)


print(symInfoFromSpacegroupAndWyckoff(10, 'd', site = None))
print(symInfoFromSpacegroupAndWyckoff(10, 'm', site = [.1, .2, .3]))


