import numpy as np
import scipy as sp
from scipy.optimize import leastsq
from os import path
from Motorcycle import GetData
gd = GetData.geda()
import matplotlib.pyplot as plt

# Inputs
#Xi = np.array([0,13,25,38,50,63,75,87,100])     # Input data X
#Yi = np.array([10.6,16.0,45.0,83.5,52.8,19.9,10.8,8.25,4.7])
#initialParameters = np.array([16000,40,400],dtype=float) # Initial Parameters of Lorentzian

file_path = path.dirname(__file__)+'/TestFile/240 50%-1.txt'
data = gd.GetRaman(file_path)
Xi = data[0]
Yi = data[1]
# InitialParameters_set = [1,640,1,2,670,2,3,690,3]
InitialParameters_set = [1,1,1,2,2,2]
# InitialParameters_set = [1,1,1]
NLorentzian = 2  # Number of Lorentzians

def Mono_Lorentzian(p,x):  # Single Lorentzian
    a1, a2, a3 = p
    return a1/((x-a2)**2+a3)

def hypothesis_func(p_set,x,n):
    f = 0
    p_list = np.copy(p_set)
    p = np.zeros((3),dtype=float)
    for i in range(n):
        p[0] = p_list[3*i]
        p[1] = p_list[3*i+1]
        p[2] = p_list[3*i+2]
        f += Mono_Lorentzian(p,x)
    return f

def Error(p_set,x,y,n):  # Measuring the Error between the hypothesis function and data
    return y-hypothesis_func(p_set,x,n)

Para = leastsq(Error,InitialParameters_set,args=(Xi,Yi,NLorentzian)) # 变量的排序需跟优化的函数中变量的排序一一对应

def R2(p_set,x,y,n):  # Coefficient of Determination
    y_average = np.sum(y)/len(y)
    SS_total = np.sum(np.square(y-y_average))
    SS_residual = np.sum(np.square(Error(p_set,x,y,n)))

    #SS_total = 0
    #SS_residual = 0
    #for i in range(len(y)):
        #SS_total += (y[i]-y_average)**2
        #SS_residual += error[i]**2

    R_squares = 1-SS_residual/SS_total
    return R_squares

def Decode(p_set,n):  # Decoding the physical properties
    property_list = np.zeros((3*n),dtype=float)
    for i in range(n):
        property_list[3*i] = p_set[3*i+1]  # Position of the maximum
        property_list[3*i+1] = 2*np.sqrt(p_set[3*i+2])  # FWHM, Full width at half maximum
        property_list[3*i+2] = p_set[3*i]/p_set[3*i+2]  # Amplitude
    return property_list

#P_test = [1,1,1]
#Para = leastsq(Error,P_test,args=(Xi,Yi))

print(Para)  # tuple
print(Para[0])
print(R2(Para[0],Xi,Yi,NLorentzian))
print(Decode(Para[0],NLorentzian))

#print(hypothesis_func(InitialParameters_set,Xi,NLorentzian)[1])
#print(Mono_Lorentzian(InitialParameters,Xi))

x0 = np.linspace(575,750,1000)
plt.scatter(Xi,Yi)
plt.plot(x0,hypothesis_func(Para[0],x0,NLorentzian),color='k')
plt.xlim(580,750)