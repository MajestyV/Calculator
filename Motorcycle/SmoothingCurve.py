import numpy as np

class smooth:
    def __init__(self):
        self.name = smooth

    def DTA_smoothing(self,x,y,n): # DTA-Discrete Trigonometric Approximation
        # pi = np.pi
        x0 = x[0]
        xn = x[len(x)]
        x_new = []
        for i in range(len(x)):
            x_new.append(x[i]/(xn-x0)-(x0+np.pi))

        a = []
        b = []
        m = len(x)/2.0
        for k in range(n):
            cos_sum = 0
            sin_sum = 0
            for j in range(len(x)):
                cos_sum += y[j]*np.cos(k*x_new[j])
                sin_sum += y[j]*np.sin(k*x_new[j])
            ak = cos_sum/m
            bk = sin_sum/m
            a.append(ak)
            if i == 0 or i == n:
                b.append(0)
            else:
                b.append(bk)

        S = []
        for i in range(len(x)):
            S0 = a[0]/2.0+a[n]*np.cos(n*x_new[i])
            sum = 0
            for k in range(1,n-1):
                sum += a[k]*np.cos(k*x[i])+b[k]*np.sin(k*x[i])
            S.append(S0+sum)

        return x,S




