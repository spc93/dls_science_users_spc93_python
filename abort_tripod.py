''' Abort script to be run if the tripod has been given commands that could damage it. '''

import platform

try:
    import argparse
except ImportError:
    if platform.system() != 'Java':
        raise

from pkg_resources import require, DistributionNotFound
try:
    require("cothread")
    from cothread.catools import caput
    from cothread.catools import caget
    from cothread.catools import DBR_CHAR_STR
except DistributionNotFound:
    if platform.system() != 'Java':
        raise

# Imports different functions depending on which computer system the code is being run on.


def main():
    caput("BL16I-MO-KBML-02:X1.STOP")
    caput("BL16I-MO-KBML-02:X2.STOP")
    caput("BL16I-MO-KBML-02:X3.STOP")
    caput("BL16I-MO-KBML-02:Y1.STOP")
    caput("BL16I-MO-KBML-02:Y2.STOP")
    caput("BL16I-MO-KBML-02:Y3.STOP")

# Stops all the motors from moving.

if __name__ == '__main__':
    main()
