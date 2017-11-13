import numpy as np

def Biaxial_Sakino_Cir(D, t, fc, fy):
    #Steel properties
    fyt = 1.08 * fy
    fyc = 0.91 * fy
    
    #Confined Concrete
    k = 4.1
    ke = 23.0
    ru = 1.67 * D**-0.112
    fcp = fc * ru
    Ec = (6.9 + 3.32 * fcp**0.5) * 1.0e3
    ec = 0.94 * fcp**0.25 * 1.0e-3
    fr = 2.0 * t * 0.19 * fy / (D - 2.0 * t)
    fre = k / ke * fr
    fcc = fcp + k * fr
    K = fcc / fcp
    if K <= 1.5:
        ecc = ec * (1.0 + 4.7 * (K - 1.0))
    else:
        ecc = ec * (3.35 + 20.0 * (K - 1.5))
    return fcc, fyc, fyt
    
    