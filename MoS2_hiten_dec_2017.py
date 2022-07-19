import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option('display.max_rows',8)

path='/dls/science/users/spc93/data/MoS2/'

file='MoS2_I16_MoS2_Point detector.txt'
d=pd.read_csv((path+file),header=250,skip_blank_lines=False)

d.plot(x='     Angle',y=' Det1Disc1', logy=True)
