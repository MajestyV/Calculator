import numpy as np
import scipy as sp
from scipy.optimize import leastsq

# Inputs
Xi = np.array([0,13,25,38,50,63,75,87,100])     # Input data X
Yi = np.array([10.6,16.0,45.0,83.5,52.8,19.9,10.8,8.25,4.7])
# InitialParameters = np.array([16000,40,400],dtype=float) # Initial Parameters of Lorentzian
InitialParameters = np.array([16000,40,400])  # leastsq can use both array or list

def func(p,x):  # Lorentzian
    a1, a2, a3 = p
    return a1/((x-a2)**2+a3)

def Error(p,x,y):  # Measuring the Error between the hypothesis function and data
    return y-func(p,x)

Para = leastsq(Error,InitialParameters,args=(Xi,Yi)) # 变量的排序需跟优化的函数中变量的排序一一对应

def R2(p,x,y):  # Coefficient of Determination
    y_average = np.sum(y)/len(y)
    error = Error(p,x,y)
    SS_total = 0
    SS_residual = 0
    for i in range(len(y)):
        SS_total += (y[i]-y_average)**2
        SS_residual += error[i]**2
    R_squares = 1-SS_residual/SS_total
    return R_squares

def Decode(p):  # Decoding the physical properties
    a1, a2, a3 = p
    FWHM = 2*np.sqrt(a3)  # Full width at half maximum
    PosMaximum = a2  # Position of the maximum
    Amplitude = a1/a3
    return PosMaximum,FWHM,Amplitude



#P_test = [1,1,1]
#Para = leastsq(Error,P_test,args=(Xi,Yi))

print(Para)
print(R2(Para[0],Xi,Yi))
print(Decode(Para[0]))