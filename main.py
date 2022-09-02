# title: Transfer Matrix Method
# Author details: Author: Furkan ÇÖZELİ, Contact details: cozelifurkan117@gmail.com
# Script and data info: This script performs transfer matrix method analyses on created multilayer structure.  
# Data consists of counts of bird species.
# Data was collected in the hunter valley region between 1990 and 1991. 
# Copyright statement: This script is the product of UNSW etc.

import numpy as np
import matplotlib.pyplot as plt
import TMMlib

class Materiel():
    """
    Since the materials have properties,
    the oop method was used.
    """
    def __init__(self, n, lambda0,  d = -1, expansionCoef = 0, opticCoef = 0, T = -1):
        """
        If the material has a custom thickness,
        the thickness value can be sent.
        Otherwise, it is calculated as follows.
        """
        self.n = n
        self.lambda0 = lambda0
        self.d = d
        if self.d == -1:
            self.d = self.lambda0 / (4 * self.n)
        if T > -1:
            self.n = n * (1 + opticCoef * (T - 300))
            self.d = d * (1 + expansionCoef * (T - 300))

mode = "TE"                             #use "TE" for TE mode, use "TM" for TM mode
lambda0 = 600 * (10**-9)
lambdaspace = np.linspace(200, 1200, 2000)*(10**-9)
angles = np.array([0])

A = Materiel(1.45, lambda0)             #SiO2
B = Materiel(4.234, lambda0)            #Te
air = Materiel(1, lambda0, 0)

"""
The thickness of the first medium and the final medium is not important.
Therefore, their thickness information needs to be entered as 0.
"""

materiels = [A, B]

series = 5*("AB")

layers = TMMlib.makeLayers(series, materiels)
structure = np.array([air])
structure = np.append(structure, layers[:])
structure = np.append(structure, air)

resultR, resultT, resultA = TMMlib.resultTMM(structure, lambdaspace, angles, mode)

#plot graph for result
plt.figure(1)

plt.subplot(3,1,1)
plt.plot(lambdaspace, resultR[0,:], label = "Reflection")
plt.ylabel("Reflection")

plt.subplot(3,1,2)
plt.plot(lambdaspace, resultT[0,:], label = "Transmission")
plt.ylabel("Transmission")

plt.subplot(3,1,3)
plt.plot(lambdaspace, resultA[0,:], label = "Absorbtion")
plt.ylabel("Absorbtion")

plt.xlabel("Wavelength")

plt.show()
