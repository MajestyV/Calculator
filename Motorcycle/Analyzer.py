import math
import numpy as np
# from numpy import *

class Lorentzian:
    def __init__(self,x_Data="",y_Data="",InitialParameters=""):
        self.x_Data = x_Data
        self.y_Data = y_Data
        self.InitialParameters = InitialParameters
        self.NumData = len(x_Data)              # Number of Data
        self.NumParm = len(InitialParameters)   # Number of Parameters
        self.name = Lorentzian

    def function(self,x,Parameters): # The Lorentzian Function
        f = Parameters[0]/((x-Parameters[1])**2+Parameters[2])
        return f


    def Error(self,a=""): # Using Least Squares Method to determine the Error function of Lorentzian funtion
        if not a:
            a = self.InitialParameters
        x = self.x_Data
        y = self.y_Data
        NumData = self.NumData
        num = self.NumParm

        Err = np.zeros((num),dtype=float) # E for Error function derived by LSM
        for i in range(NumData):
            Err[0] += (self.function(x[i],a)-y[i])/((x[i]-a[1])**2+a[2])
            Err[1] += (self.function(x[i],a)-y[i])*(-2*a[0]*(x[i]-a[2]))/((x[i]-a[1])**2+a[2])**2
            Err[2] += (self.function(x[i],a)-y[i])*(-a[0])/((x[i]-a[1])**2+a[2])**2

        return Err

    def Jacobian(self,a="",da="",num=""): # Determine the Jacobian used in Newton's Method
        if not a:
            a = self.InitialParameters
        if not da:
            da = 0.00001
        #if not num:
            #num = self.NumData
        J = np.zeros((num,num),dtype=float)
        a_plus = np.copy(a)
        a_minus = np.copy(a)

        for i in range(num):
            for j in range(num):
                a_plus[j] = a[j]+da/2      # x+dx/2
                a_minus[j] = a[j]-da/2     # x-dx/2
                J[i,j] = (self.Error(a_plus)[i]-self.Error(a_minus)[i])/da # [f(x+dx/2)-f(x-dx/2)]/dx

        return J

    def Newton(self,a=""): # Using Newton's Method to optimize parameters
        if not a:
            a = self.InitialParameters
        a1 = np.copy(a)
        InvJ = np.linalg.inv(self.Jacobian(a)) # Inverse of Jacobian
        Error = np.copy(a)
        InvJ0 = InvJ
        i = 0
        while(np.sum(abs(Error))>1.e-6 and i<100):
            a1 = a-np.dot(np.linalg.inv(self.Jacobian(a)))



