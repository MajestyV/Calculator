import numpy as np
import matplotlib.pyplot as plt

N = 20
e = 1.6*10**-19
m = 9.109*10**-31
c = 3*10**8
R = 0.25
I0 = 100000
N1 = 1000
N2 = 1000
N3 = 1000

def XRD(theta,CrystalPlane):
    y1 = I0*e**4*(1+np.cos(2*theta)**2)/2*m**2*c**4*R**2
    y2 = func(CrystalPlane)
    return y1*y2

def func(CrystalPlane):
    h, k, l = CrystalPlane
    if h == 0: h = 0.00001
    if k == 0: k = 0.00001
    if l == 0: l = 0.00001
    f = np.sin(N1 * np.pi * h) ** 2 * np.sin(N2 * np.pi * k) ** 2 * np.sin(N3 * np.pi * l) ** 2 / (np.sin(np.pi * h) ** 2 * np.sin(np.pi * k) ** 2 * np.sin(np.pi * l) ** 2)
    return f

#x = np.linspace(0.0001,360,1000)
#plt.plot(x,XRD(x,[1,1,8]))
#plt.ylim(-1,10)
#plt.show()

#print(func([0,0,1]))

#print(np.array([4,6,8])*np.array([7,23,4]))
#print([4,6,1]-[1,2,4])

def StructureFactor(AtomicPosition,AtomicScatteringFactor,CrystalPlane):
    FF_star_list = []
    for n in range(len(CrystalPlane)):
        g = np.array(CrystalPlane[n])
        ASF = [AtomicScatteringFactor[m][n] for m in range(len(AtomicScatteringFactor))]  # ASF - Atomic Scattering Factor
        FF_star = 0
        for i in range(len(AtomicPosition)):
            for j in range(len(AtomicPosition)):
                ra = np.array(AtomicPosition[i])  # Atomic Position Vector
                rb = np.array(AtomicPosition[j])
                fa = ASF[i]  # Atomic Scattering Factor
                fb = ASF[j]
                FF_star += fa*fb*np.exp(2*np.pi*(0.0+1j)*np.sum(g*(ra-rb)))
        FF_star_list.append(FF_star.real)
    return FF_star_list

# [[np.array([0,0,0]),1], [np.array([0.5,0.5,0.5]),1]]
AP = [[0,0,0],[0.5,0.5,0.5]]
ASF = [[1,1,1,1,1,1],[1,1,1,1,1,1]]
HKL = [[1,0,0],[2,0,0],[3,0,0],[4,0,0],[5,0,0],[6,0,0]]

print(StructureFactor(AP,ASF,HKL))