import numpy as np
import pycbc.inject


def readMasses(file):
    injections = pycbc.inject.InjectionSet(file)
    mass1 = []
    mass2 = []
    for ii, inj in enumerate(injections.table):
        mass1.append(inj.mass1)
        mass2.append(inj.mass2)
    mass1 = np.array(mass1)
    mass2 = np.array(mass2)
    masses = np.vstack((mass1, mass2)).T
    return masses
    
    
    