import numpy as np
import sys



def rgwlc(Lp,Lc):
    rg2 = (1.0/3.0)*Lp*Lc - Lp**2 + 2*((Lp**3) / Lc)* (1- (Lp/Lc)*(1- np.exp(-Lc/Lp)))

    return np.sqrt(rg2)


Lp = 10.0
print rgwlc(Lp,float(sys.argv[1]))