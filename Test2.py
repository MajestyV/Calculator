import math
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import matplotlib

#plt.figure(figsize=(20, 16), dpi=100)
#plt.rc('font', family='Times New Roman')
#matplotlib.rcParams['font.family'] = 'STSong'
#plt.rcParams['axes.unicode_minus'] = False
#matplotlib.rcParams['font.size'] = 20

ND = 9
num = 3
E = np.array([0, 13, 25, 38, 50, 63, 75, 87, 100])
S = np.array([10.6, 16.0, 45.0, 83.5, 52.8, 19.9, 10.8, 8.25, 4.7])
sig = np.ones_like(S)
# sig =np.array([9.34, 17.9, 41.5, 85.5, 51.5, 21.5, 10.8, 6.29, 4.14])    # Error bars
x = np.array([20000, 40, 220], dtype=float)
#plt.plot(E, S, 'bo', markersize=20, alpha=0.4, label="采样数据点")

#plt.title('Least-Squares Fit of linear to Blue Data', fontsize=20)
#plt.xlabel('E', fontsize=32);
#plt.ylabel('$\sigma$', fontsize=32);
#plt.grid(True)


def Fun(x, num):
    F = np.zeros((num), dtype=float)
    for j in range(ND):
        F[0] += (S[j] - x[0] / ((E[j] - x[1]) ** 2 + x[2])) / sig[j] ** 2 / ((E[j] - x[1]) ** 2 + x[2])
        F[1] += -2 * x[0] * (E[j] - x[1]) * (S[j] - x[0] / ((E[j] - x[1]) ** 2 + x[2])) / sig[j] ** 2 / ((E[j] - x[1]) ** 2 + x[2]) ** 2
        F[2] += -x[0] * (S[j] - x[0] / ((E[j] - x[1]) ** 2 + x[2])) / sig[j] ** 2 / ((E[j] - x[1]) ** 2 + x[2]) ** 2

    return F


def dfun(x, num):
    df = np.zeros((num, num), dtype=float)
    dx = 0.00001  #
    x1_plus = np.copy(x)
    x1_minus = np.copy(x)

    for i in range(0, num):
        for j in range(0, num):
            x1_plus = np.copy(x)
            x1_minus = np.copy(x)
            x1_plus[j] = x1_plus[j] + dx / 2  # x+dx
            x1_minus[j] = x1_minus[j] - dx / 2
            df[i, j] = (Fun(x1_plus, num)[i] - Fun(x1_minus, num)[i]) / dx  # f(x+dx)-f（x）/dx
    df_1 = np.linalg.inv(df)  # 计算逆矩阵
    return df_1


def Newton(x, num):
    x1 = np.copy(x)
    i = 0
    delta = np.copy(x)
    #    dfun0=dfun(x,num)
    while (np.sum(abs(delta)) > 1.e-6 and i < 100):  # 控制循环次数
        x1 = x - np.dot(dfun(x, num), Fun(x, num))  # 公式
        delta = x1 - x
        x = x1
        i = i + 1
        print(x, )
    return x


a = Newton(x, num)
print(a)

#EX = np.linspace(0, 100, 1000)
#G_E = x[0] / ((EX - x[1]) ** 2 + x[2])
#G = a[0] / ((EX - a[1]) ** 2 + a[2])
#points = a[0] / ((E - a[1]) ** 2 + a[2])
#plt.plot(EX, G_E, 'k', linewidth=3.0, label="初猜曲线")
#plt.plot(EX, G, 'r:', linewidth=5.0, label="最终拟合曲线")
#plt.plot(E, points, 'ro', markersize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#plt.legend()
#plt.show()