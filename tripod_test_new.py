from tripod_test_class import *

from pkg_resources import require, DistributionNotFound
try:
    require("cothread")
    from cothread.catools import caput
    from cothread.catools import caget
    from cothread.catools import DBR_CHAR_STR
except DistributionNotFound:
    if platform.system() != 'Java':
        raise

tp = tripod_class([134.2, 134.2, 134.2], [219.129, 219.129, 84.963], [-pi/3, pi/3, 0], [150.102, 84.9634/2, 35.7574], [pi/4, pi/4, -pi/4], [0.0, 0.0, 357.313], [249.324, 0.0, 249.324/2])

tripod=tp


X, Y = [0, 0, 0], [0, 0, 0]
C, alpha = tripod.ctool(X, Y)


class stop_all():
	def __init__(self):
		print '===  IMPORTANT: Hit q,a or z <return> to abort motors  ==='
	def __repr__(self):
		caput("BL16I-MO-KBML-02:X1.STOP",1)
		caput("BL16I-MO-KBML-02:X2.STOP",1)
		caput("BL16I-MO-KBML-02:X3.STOP",1)
		caput("BL16I-MO-KBML-02:Y1.STOP",1)
		caput("BL16I-MO-KBML-02:Y2.STOP",1)
		caput("BL16I-MO-KBML-02:Y3.STOP",1)
		return '=== All motors asked to stop'
a=q=z=stop_all()

Czero=np.array([0,0,0])
alphazero=np.array([0,0,0])

def ctool():
	'Returns current tripod settings'
	X[0]=caget("BL16I-MO-KBML-02:X1.RBV")
	X[1]=caget("BL16I-MO-KBML-02:X2.RBV")
	X[2]=caget("BL16I-MO-KBML-02:X3.RBV")
	Y[0]=caget("BL16I-MO-KBML-02:Y1.RBV")
	Y[1]=caget("BL16I-MO-KBML-02:Y2.RBV")
	Y[2]=caget("BL16I-MO-KBML-02:Y3.RBV")
	C, alpha = tripod.ctool(X, Y)
	return (C, alpha)

def show():
	ctool()
	return tp.__repr__
	

def mv_angs(alpha_new):
	'Move tripod to new angles keeping translations the same'
	C, alpha = ctool()
	X, Y = tripod.cbase(C, np.array(alpha_new))
	caput("BL16I-MO-KBML-02:X1.VAL", (X[0]))
	caput("BL16I-MO-KBML-02:X2.VAL", (X[1]))
	caput("BL16I-MO-KBML-02:X3.VAL", (X[2]))
	caput("BL16I-MO-KBML-02:Y1.VAL", (Y[0]))
	caput("BL16I-MO-KBML-02:Y2.VAL", (Y[1]))
	caput("BL16I-MO-KBML-02:Y3.VAL", Y[2])

def mv_trans(Cnew):
	'Move tripod to new translations keeping angles the same'
	C, alpha = ctool()
	X, Y = tripod.cbase(Cnew, alpha)
	caput("BL16I-MO-KBML-02:X1.VAL", (X[0]))
	caput("BL16I-MO-KBML-02:X2.VAL", (X[1]))
	caput("BL16I-MO-KBML-02:X3.VAL", (X[2]))
	caput("BL16I-MO-KBML-02:Y1.VAL", (Y[0]))
	caput("BL16I-MO-KBML-02:Y2.VAL", (Y[1]))
	caput("BL16I-MO-KBML-02:Y3.VAL", Y[2])









