import numpy as np
import matplotlib.pyplot as plt

# 固有常数及机器参数
e = 1.6*10**-19
m = 9.109*10**-31
c = 3*10**8
R = 0.25
I0 = 100000

# 样品参数
N1 = 10000
N2 = 10000
N3 = 10000
SampleInformation = [[np.array([0,0,0]),1], [np.array([0.5,0.5,0.5]),1]]

#测试参数
def CrystalPlaneFamily(h_max,k_max,l_max):
    hkl = []
    for h in range(h_max+1):
        for k in range(k_max+1):
            for l in range(l_max+1):
                if [h,k,l] != [0,0,0]:
                    hkl.append([h,k,l])
    return hkl

def Intensity(Theta):
    Ie = I0*e**4*(1+np.cos(2*Theta)**2)/2*m**2*c**4*R**2
    return Ie

def SampleFactor(CrystalPlane):
    sf_list = []  # abbreviation for sample factor
    for i in range(len(CrystalPlane)):
        h, k, l = CrystalPlane[i]
        if h == 0: h = 0.00001
        if k == 0: k = 0.00001
        if l == 0: l = 0.00001
        f = np.sin(N1*np.pi*h)**2*np.sin(N2*np.pi*k)**2*np.sin(N3*np.pi*l)**2/(np.sin(np.pi*h)**2*np.sin(np.pi*k)**2*np.sin(np.pi*l)**2)
        sf_list.append(f)
    return sf_list

def StructureFactor(StructureInformation,CrystalPlane):
    FF_star_list = []
    for n in range(len(CrystalPlane)):
        g = np.array(CrystalPlane[n])
        FF_star = 0
        for i in range(len(StructureInformation)):
            for j in range(len(StructureInformation)):
                ra = StructureInformation[i][0]  # Atomic Position Vector
                rb = StructureInformation[j][0]
                fa = StructureInformation[i][1]  # Atomic Scattering Factor
                fb = StructureInformation[j][1]
                FF_star += fa*fb*np.exp(2*np.pi*(0.0+1j)*np.sum(g*(ra-rb)))
        FF_star_list.append(FF_star.real)
    return FF_star_list

# 主函数
def XRD(Theta,StructureInformation,h_max,k_max,l_max):
    CrystalPlane = CrystalPlaneFamily(h_max,k_max,l_max)
    SamFac = SampleFactor(CrystalPlane)
    StrucFac = StructureFactor(StructureInformation,CrystalPlane)

    I_list = []
    Isum = []
    Isum_list = []
    for i in range(len(Theta)):
        for j in range(len(CrystalPlane)):
            I = Intensity(Theta[i])*SamFac[j]*StrucFac[j]
            I_list.append(I)
        Isum = np.sum(I_list)
        Isum_list.append(Isum)

    return Isum_list







#def Intensity(Theta,CrystalPlane,N1=10000,N2=10000,N3=10000):
    #Ie = I0*e**4*(1+np.cos(2*Theta)**2)/2*m**2*c**4*R**2

a = CrystalPlaneFamily(5,5,5)
print(a)
print(len(a))

print(StructureFactor(SampleInformation,CrystalPlaneFamily(0,0,2)))

print(np.sum([1,1,1]))
print(XRD([20,40,60],SampleInformation,3,3,3))
#print(len(XRD(20,SampleInformation,3,3,3)))

#print(SampleFactor([[0,0,1]]))

doublet = np.linspace(0,90,1000)

plt.plot(doublet,XRD(doublet/2,SampleInformation,3,3,3))
plt.ylim(-1,10)
plt.show()