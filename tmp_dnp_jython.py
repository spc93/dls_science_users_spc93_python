
#import scisoftpy as np
import numpy as np


for i in range(100):
    print i
    for i in range(10000):    
        a6=np.random.rand(6,6)
        a6i=np.linalg.inv(a6)

np.dot(a6i, a6)

# does approx. 10K inversions/sec - not much faster in numpy