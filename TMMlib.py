import matplotlib.pyplot as plt
import numpy as np

def makeLayers(series, materiels):
    layers = np.array([materiels[ord(series[0])-65]])
    series = series[1:]

    for layer in series:
        layers = np.append(layers, materiels[ord(layer)-65])
        
    return layers

def transferMat(rho, tao, phi):
    # Transfer Matrix 2x2 matrix
    return np.array([[np.exp(1j*phi), rho*np.exp(-1j*phi)],
                     [rho*np.exp(1j*phi), np.exp(-1j*phi)]])/tao

def taoTM(n1, n2, teta1, teta2):
    #parallel polarization
    return (2*n1*np.cos(teta1)) / (n1*np.cos(teta2) + n2*np.cos(teta1))

def taoTE(n1, n2, teta1, teta2):
    #perpendicular polarization
    return (2*n1*np.cos(teta1)) / (n1*np.cos(teta1) + n2*np.cos(teta2))

def rhoTM(n1, n2, teta1, teta2):
    #parallel polarization
    return (n1*np.cos(teta2)-n2*np.cos(teta1)) / (n1*np.cos(teta2)+n2*np.cos(teta1))

def rhoTE(n1, n2, teta1, teta2):
    #perpendicular polarization
    return (n1*np.cos(teta1)-n2*np.cos(teta2)) / (n1*np.cos(teta1)+n2*np.cos(teta2))

def teta(n1, n2, teta1):
    return np.arcsin(n1 * np.sin(teta1) / n2)

def phi(wavelength, d, n, teta):
    return 2 * np.pi * d * n * np.cos(teta) / wavelength

def TMM(mode, wavelength, structure, angle):
    """
    Calculate the behavior of the multilayer structure for the respective angle and wavelength.
    Both TE and TM polarizations have been added separately.
    """
    sumOfLayers = 0
    teta1 = np.radians(angle)
    if mode == "TE":
        i = 0
        while i < structure.size - 1:
            n1, n2 = structure[i].n, structure[i+1].n
            teta2 = teta(n1, n2, teta1)
            d = structure[i+1].d
            rho = rhoTE(n1, n2, teta1, teta2)
            tao = taoTE(n1, n2, teta1, teta2)
            calPhi = phi(wavelength, d, n2, teta2)
            teta1 = teta2
            
            if i == 0:
                transferMatrix = transferMat(rho, tao, calPhi)
                i += 1
                continue
            transferMatrix = np.matmul(transferMatrix, transferMat(rho, tao, calPhi))
            i += 1
        
        R = np.absolute(transferMatrix[1, 0] / transferMatrix[0, 0])**2
        T = (structure[structure.size-1].n / structure[0].n) * np.absolute(1 / transferMatrix[0, 0])**2
        A = 1 - R - T
        return R, T, A
        

    if mode == "TM":
        i = 0
        while i < structure.size -1:
            n1, n2 = structure[i].n, structure[i+1].n
            teta2 = teta(n1, n2, teta1)
            d = structure[i+1].d
            rho = rhoTM(n1, n2, teta1, teta2)
            tao = taoTM(n1, n2, teta1, teta2)
            calPhi = phi(wavelength, d, n2, teta2)
            teta1 = teta2
            
            if i == 0:
                transferMatrix = transferMat(rho, tao, calPhi)
                i += 1
                continue
            transferMatrix = np.matmul(transferMatrix, transferMat(rho, tao, calPhi))
            i += 1
        
        R = np.absolute(transferMatrix[1, 0] / transferMatrix[0, 0])**2
        T = (structure[structure.size-1].n / structure[0].n) * np.absolute(1 / transferMatrix[0, 0])**2
        A = 1 - R - T
        return R, T, A

def resultTMM(structure, lambdaspace, angles, mode):
    """
    Here, the final result is obtained for both multiple angles and wavelength range.
    Each line calculates the multilayer structure behavior in the specific wavelength
    range for the respective angle.
    """
    resultR = np.zeros((angles.size, lambdaspace.size))
    resultT = np.zeros((angles.size, lambdaspace.size))
    resultA = np.zeros((angles.size, lambdaspace.size))
    i = 0

    while i < angles.size:
        j = 0

        while j < lambdaspace.size:
            reflectivity, transmissivity, absorption = TMM(mode,lambdaspace[j], structure, angles[i])

            if j == 0:
                reflectivityMat = np.array([reflectivity])
                transmissivityMat = np.array([transmissivity])
                absorptionMat = np.array([absorption])
                j += 1
                continue

            reflectivityMat = np.append(reflectivityMat, reflectivity)
            transmissivityMat = np.append(transmissivityMat, transmissivity)
            absorptionMat = np.append(absorptionMat, absorption)
            j += 1

        resultR[i, : ] = reflectivityMat[:]
        resultT[i, : ] = transmissivityMat[:]
        resultA[i, : ] = absorptionMat[:]
        i += 1

    return resultR, resultT, resultA
