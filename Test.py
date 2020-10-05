import numpy as np
from scipy.optimize import leastsq
import pylab as pl
import matplotlib.pyplot as plt

def func(x, p):
    a1, a2, a3 = p
    return a1/((x-a2)**2+a3)   

def residuals(p, y, x):
    return y - func(x, p)

x = np.linspace(0, 100, 1000)

p0 = [16000,40, 200, ]
x0=np.array([0, 13 ,25 ,38, 50, 63, 75, 87, 100])
y0=np.array([10.6, 16.0, 45.0, 83.5, 52.8, 19.9, 10.8, 8.25, 4.7])
y1 = y0 + 2 * np.random.randn(len(x0))

plsq = leastsq(residuals, p0, args=(y1, x0))

print (u"拟合参数", plsq[0])
plt.title(u'洛伦兹函数最小二乘法拟合',fontsize=20)
pl.plot(x0, y0, 'bo',markersize=20,alpha=0.4, label=u"真实数据")
pl.plot(x0, y1, 'y*', markersize=20,alpha=0.4, label=u"带噪声的实验数据")
pl.plot(x, func(x, plsq[0]), 'r:',linewidth=5.0, label=u"拟合曲线")
plt.plot(x, func(x, p0),'k', linewidth=3.0,label="初猜曲线")
plt.plot(x0, func(x0, plsq[0]), 'ro',markersize= 20)
plt.xlabel('E',fontsize=32);  plt.ylabel('$\sigma$',fontsize=32);   plt.grid(True)
pl.legend()
pl.show()