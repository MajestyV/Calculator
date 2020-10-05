import numpy as np
import matplotlib.pyplot as plt

# 输入参数测试
# [[Atomic Coordinate], Atomic Number/Element]
StrucInform = [[[0,0,0],17],[[0.5,0.5,0],17],[[0,0.5,0.5],17],[[0.5,0,0.5],17],
               [[0.5, 0.5, 0.5],11],[[0,0,0.5],11],[[0,0.5,0],11],[[0.5,0,0],11]]
#StrucInform = [[[0,0,0],28], [[0.5,0.5,0],28], [[0.5,0,0.5],28], [[0,0.5,0.5],28]]
#LP = [3.5288,3.5288,3.5288,90,90,90,'Cubic']
LP = [5.6071,5.6071,5.6071,90,90,90,'Cubic']


# 0-AtomName, 1-AtomicNumber, 2-a1, 3-b1, 4-a2, 5-b2, 6-a3, 7-b3, 8-a4, 9-b4
AtomicScatteringInformation = [['Ac', 98, 6.278, 28.323, 5.195, 4.949, 2.321, 0.557, 0, 0],
                               ['Ag', 47, 2.036, 61.497, 3.272, 11.824, 2.511, 2.846, 0.837, 0.327],
                               ['Na', 11, 2.241, 108.004, 1.333, 24.505, 0.907, 3.391, 0.286, 0.435],
                               ['Cl', 17, 1.452, 30.935, 2.292, 9.980, 0.787, 2.234, 0.322, 0.323],
                               ['W', 74, 5.709, 28.782, 4.677, 5.084, 2.019, 0.572, 0, 0],
                               ['Ni',28, 2.210, 58.727, 2.134, 13,553, 1.689, 2.609, 0.524, 0.339]]

# 波形拟合函数组
def KroneckerDelta(i,j):
    y = 0
    for n in range(len(j)):
        if i == j[n]:
            y = 1
    return y


def DiracDelta(i,j,a=0.00001):
    return np.exp(-(i-j)**2/a**2)/(a*np.sqrt(np.pi))

# 峰位计算函数组
def ReciprocalMetricTensor(a,b,c,alpha,beta,gamma,CrystalSystem='Unknown'):
    if CrystalSystem == 'Unknown':
        a_star = a
        b_star = b
        c_star = c
        alpha_star = alpha
        beta_star = beta
        gamma_star = gamma
        g_star = np.array([[a_star**2, a_star*b_star*np.cos(gamma_star), a_star*c_star*np.cos(beta_star)],
                           [b_star*a_star*np.cos(gamma_star), b_star**2, b_star*c_star*np.cos(alpha_star)],
                           [c_star*a_star*np.cos(beta_star), c_star*b_star*np.cos(alpha_star), c_star**2]])
    elif CrystalSystem == 'Cubic':
        g_star = np.array([[1/a**2, 0, 0],
                           [0, 1/a**2, 0],
                           [0, 0, 1/a**2]])
    elif CrystalSystem == 'Orthorhombic':
        g_star = np.array([[1/a**2, 0, 0],
                           [0, 1/b**2, 0],
                           [0, 0, 1/c**2]])

    return g_star

def CrystalPlaneFamily(h_max,k_max,l_max):
    hkl = []
    for h in range(h_max+1):
        for k in range(k_max+1):
            for l in range(l_max+1):
                if [h,k,l] != [0,0,0]:
                    hkl.append([h,k,l])
    return hkl

def PlaneToAngle(CrystalPlane,Wavelength,LatticeParameters):  # Converting a set of Crystal Plane Families with Miller Index to a set of angles(theta).
    a, b, c, alpha, beta, gamma, CrystalSystem = LatticeParameters
    g_star = ReciprocalMetricTensor(a,b,c,alpha,beta,gamma, CrystalSystem)

    hkl_list = []
    dhkl_list = []
    OutOfRange_list = []
    Theta_list = []
    for i in range(len(CrystalPlane)):
        hkl = np.array(CrystalPlane[i])
        dhkl = 1/np.sqrt(np.sum(hkl*g_star*np.transpose(hkl)))
        if Wavelength/(2*dhkl) <= 1:
            Theta = np.arcsin(Wavelength/(2*dhkl))*180/np.pi
            Theta_list.append(Theta)
            hkl_list.append(hkl)
            dhkl_list.append(dhkl)
        else:
            OutOfRange_list.append(hkl)

    if OutOfRange_list == []:
        print('All the Crystal Plane Families are inside the test range.')
    else:
        print('The below Crystal Plane Families are out of range:')
        print(OutOfRange_list)

    return Theta_list,hkl_list,dhkl_list

def AtomicScatteringFactor(Theta,Wavelength,Z,T=0,AtomicScatteringParameters=AtomicScatteringInformation):
    s = np.sin(Theta)/Wavelength
    a = np.zeros((4),dtype=float)
    b = np.zeros((4),dtype=float)
    sum = 0
    for i in range(len(AtomicScatteringParameters)):
        if Z == AtomicScatteringParameters[i][1]:
            for n in range(4):
                a[n] = AtomicScatteringParameters[i][2*n+2]
                b[n] = AtomicScatteringParameters[i][2*n+3]
                sum += a[n]*np.exp(-b[n]*s**2)

    # Debye-Waller Temperature Factor
    #TemFac =

    f = Z-42.78214*s**2*sum
    return f

def StructureFactor(AtomicPosition, AtomScatFac, CrystalPlane):
    FF_star_list = []
    for n in range(len(CrystalPlane)):
        g = np.array(CrystalPlane[n])
        ASF = [AtomScatFac[m][n] for m in range(len(AtomScatFac))]  # ASF - Atomic Scattering Factor
        FF_star = 0
        for i in range(len(AtomicPosition)):
            for j in range(len(AtomicPosition)):
                ra = np.array(AtomicPosition[i])  # Atomic Position Vector
                rb = np.array(AtomicPosition[j])
                fa = ASF[i]  # Atomic Scattering Factor
                fb = ASF[j]
                FF_star += fa * fb * np.exp(2 * np.pi * (0.0 + 1j) * np.sum(g * (ra - rb)))
        FF_star_list.append(FF_star.real)
    return FF_star_list

# 主函数
def Test():
    wavelength = 1.789
    Theta_list = []
    dhkl_list = []
    hkl = CrystalPlaneFamily(5,5,5)
    hkl_list = []
    # FF_star = StructureFactor(StructureInformation,hkl)
    for i in range(len(hkl)):
        dhkl = 1/np.sqrt(np.sum(np.array(hkl[i])*ReciprocalMetricTensor(5.6071,5.6071,5.6071,90,90,90,'Cubic')*np.transpose(np.array(hkl[i]))))  # 需改进！！！
        dhkl_list.append(dhkl)
        Theta = np.arcsin(wavelength/(2*dhkl))
        Theta_list.append(Theta*180/np.pi)
        hkl_list.append(hkl[i])
    return dhkl_list,Theta_list,hkl_list

def main(StructureInformation,LaParam,TestRange):
    wavelength = 1.789

    hmax, kmax, lmax = TestRange
    hkl_old = CrystalPlaneFamily(hmax,kmax,lmax)
    theta, hkl_new, dhkl = PlaneToAngle(hkl_old,wavelength,LaParam)
    # double_theta = [theta[n]*2 for n in range(len(theta))]

    AtomPos = [StructureInformation[n][0] for n in range(len(StructureInformation))]  # Atomic Position list
    AtomNum = [StructureInformation[n][1] for n in range(len(StructureInformation))]  # Atomic Number list
    AtomScatFac = [AtomicScatteringFactor(theta,wavelength,AtomNum[n]) for n in range(len(AtomNum))]  # Atomic Scattering Factor list

    StrucFac = StructureFactor(AtomPos,AtomScatFac,hkl_new)

    theta_list = []
    # double_theta_list = []
    hkl_list = []
    dhkl_list = []
    StrucFac_new = []
    for i in range(len(StrucFac)):
        if StrucFac[i] != 0:
            if StrucFac[i] not in StrucFac_new:
                StrucFac_new.append((StrucFac[i]))
                theta_list.append(theta[i])
                # double_theta_list.append((double_theta[i]))
                hkl_list.append((hkl_new[i]))
                dhkl_list.append(dhkl[i])
            else:
                a = StrucFac_new.index(StrucFac[i])
                if not isinstance(hkl_list[a],list):
                    hkl_list[a] = [hkl_list[a]]
                hkl_list[a].append(hkl_new[i])

    theta_rearrange = sorted(theta_list)
    index_list = [theta_list.index(theta_rearrange[n]) for n in range(len(theta_rearrange))]
    hkl_rearrange = [hkl_list[index_list[n]] for n in range(len(index_list))]
    StrucFac_rearrange = [StrucFac_new[index_list[n]] for n in range(len(index_list))]
    dhkl_rearrange = [dhkl_list[index_list[n]] for n in range(len(index_list))]
    double_theta_list = [2*theta_rearrange[n] for n in range(len(theta_rearrange))]

    #print(StrucFac_new)
    #print(theta_list)
    #print(double_theta_list)
    #print(hkl_list)
    #print(theta_rearrange)
    #print(hkl_rearrange)

    return theta_rearrange,double_theta_list,StrucFac_rearrange,hkl_rearrange,dhkl_rearrange

# 可视化模块
def XRD(StartingPoint,StoppingPoint,step,double_theta,intensity):
    n = len(double_theta)
    x0 = np.linspace(0,double_theta[0],int(double_theta[0]/step),endpoint=False)
    for i in range(n-1):
        x0 = np.hstack((x0,np.linspace(double_theta[i],double_theta[i+1],int((double_theta[i+1]-double_theta[i])/step),endpoint=False)))
    x0 = np.hstack((x0,np.linspace(double_theta[n-1],360,int((360-double_theta[n-1])/step),endpoint=False)))

    x = []
    for i in range(len(x0)):
        if x0[i] >= StartingPoint and x0[i] <= StoppingPoint:
            x.append(x0[i])

    #y = []
    #for i in range(len(x)):
        #k1,k2 = KroneckerDelta(x[i],double_theta)
        #y.append(k1*intensity(double_theta.index(k2)))

    y =[]
    for i in range(len(x)):
        y0 = 0
        I = 0
        for j in range(len(double_theta)):
            if x[i] == double_theta[j]:
                y0 = 1
                I = intensity[j]
        y.append(y0*I)

    print(x)
    print(y)

    return x,y




#x = np.hstack((np.linspace(-10,0,1000),np.linspace(0.00001,10,1000)))
#plt.plot(x,KroneckerDelta(x,0),color='k')
#plt.plot(x,DiracDelta(x,0),color='y')
#plt.show()

#print(Test()[0])
#print(Test()[1])
#print(Test()[2])
#print(Test()[3])
#print(np.arcsin(0.7)*180/np.pi)
#print(np.sqrt(np.sum(np.array([1,0,0])*ReciprocalMetricTensor(5,5,5,90,90,90,'Cubic')*np.transpose(np.array([1,0,0])))))
#print(np.sqrt(np.sum(np.transpose(np.array([1,0,0]))*ReciprocalMetricTensor(5,5,5,90,90,90,'Cubic')*np.array([1,0,0]))))

#print(f(20,1.789,11))
#print(f(x,1.780,11))


#print([StrucInform[n][0] for n in range(len(StrucInform))])
a = main(StrucInform,LP,[5,5,5])

xn,yn = XRD(20,90,0.05,a[1],a[2])

plt.plot(xn,yn)
# plt.ylim(-10,2000)
for i in range(len(a[1])):
    if isinstance(a[3][i],list):
        plt.text(a[1][i]-3,a[2][i]+200,str(a[3][i][1]))
    else:
        plt.text(a[1][i]-3, a[2][i] + 200, str(a[3][i]))
