import numpy as np

# Inputs
Xi = np.array([0,13,25,38,50,63,75,87,100])     # Input data X
Yi = np.array([10.6,16.0,45.0,83.5,52.8,19.9,10.8,8.25,4.7])
sig = np.ones_like(Yi)
# InitialParameters = np.array([16000,40,400],dtype=float) # Initial Parameters of Lorentzian
InitialParameters = [20000,40,220]  # 对于牛顿法而言，初始值（即迭代起点）是否关键，不然很容易不收敛
# InitialParameters = [1,1,1]

#dp = 0.00001  # 计算Jacobian时的自变量步长

Tolerance = 1.e-6  # 牛顿法的迭代过程控制量
IterationMaximum = 1000

def func(p,x):
    a1, a2, a3 = p
    return a1/((x-a2)**2+a3)

#def func(p,x):
    #a1 = p[0]
    #a2 = p[1]
    #a3 = p[2]
    #return a1/((x-a2)**2+a3)

def ExtremumCondition(p,x,y):  # 极值条件
    a1, a2, a3 = p
    #a1 = p[0]
    #a2 = p[1]
    #a3 = p[2]
    EC = np.zeros((len(p)),dtype=float)  # Abbreviation for Extremum Condition
    for i in range(len(x)):
        # EC[0] += (func(x[i],p)-y[i])/((x[i]-a2)**2+a3)  # [h(Xi)-Yi]*dh(Xi)/da0
        # EC[1] += (func(x[i],p)-y[i])*(-2*a1*(x[i]-a3))/((x[i]-a2)**2+a3)**2
        # EC[2] += (func(x[i],p)-y[i])*(-a1)/((x[i]-a2)**2+a3)**2

        #EC[0] += (func(p,x[i]) - y[i]) / ((x[i] - a2) ** 2 + a3)  # [h(Xi)-Yi]*dh(Xi)/da0
        #EC[1] += (func(p,x[i]) - y[i]) * (-2 * a1 * (x[i] - a3)) / ((x[i] - a2) ** 2 + a3) ** 2
        #EC[2] += (func(p,x[i]) - y[i]) * (-a1) / ((x[i] - a2) ** 2 + a3) ** 2

        EC[0] += (y[i]-a1/((x[i] - a2) ** 2 + a3)) / ((x[i] - a2) ** 2 + a3)  # [h(Xi)-Yi]*dh(Xi)/da0
        EC[1] += (y[i]-a1/((x[i] - a2) ** 2 + a3)) * (-2 * a1 * (x[i] - a2)) / ((x[i] - a2) ** 2 + a3) ** 2
        EC[2] += (y[i]-a1/((x[i] - a2) ** 2 + a3)) * (-a1) / ((x[i] - a2) ** 2 + a3) ** 2

    return EC

def Fun(p,x,y):
    a1, a2, a3 = p
    F = np.zeros((len(p)), dtype=float)
    for j in range(len(x)):
        F[0] += (y[j] - a1 / ((x[j] - a2) ** 2 + a3)) / sig[j] ** 2 / ((x[j] - a2) ** 2 + a3)
        F[1] += -2 * a1 * (x[j] - a2) * (y[j] - a1 / ((x[j] - a2) ** 2 + a3)) / sig[j] ** 2 / ((x[j] - a2) ** 2 + a3) ** 2
        F[2] += -a1 * (y[j] - a1 / ((x[j] - a2) ** 2 + a3)) / sig[j] ** 2 / ((x[j] - a2) ** 2 + a3) ** 2

    return F

def Jacobian(p,x,y):
    dp =0.00001
    p_plus = np.copy(p)
    p_minus = np.copy(p)
    # a = ExtremumCondition(p,x,y)
    J = np.zeros((len(p),len(p)),dtype=float)
    for i in range(len(p)):
        for j in range(len(p)):
            #p_plus = p
            #P_minus = p
            p_plus[j] = p_plus[j]+dp/2
            p_minus[j] = p_minus[j]-dp/2
            J[i,j] = (ExtremumCondition(p_plus,x,y)[i]-ExtremumCondition(p_minus,x,y)[i])/dp

            # J[i, j] = (Fun(p_plus,x,y)[i] - Fun(p_minus,x,y)[i]) / dp  # f(x+dx)-f（x）/dx

            p_plus = np.copy(p)  # 复位
            p_minus = np.copy(p)

    # InvJ = np.linalg.inv(J)  # Inverse of Jacobian

    return J

def Newton(p,x,y):
    p_list = []
    p_list.append(p)

    p1 = np.copy(p)
    p2 = np.copy(p)
    error = np.ones_like(p)
    i = 0
    while(np.sum(np.abs(error)) > Tolerance and i < IterationMaximum):
        InJ = np.linalg.inv(Jacobian(p1,x,y))  # Inverse of Jacobian
        p2 = p1-np.dot(InJ,Fun(p1,x,y))
        error = p2-p1
        p1 = p2
        i = i+1
        print(p2)
        p_list.append(p2)

    #p2 = p1-np.dot(Jacobian(p1, x, y), ExtremumCondition(p1, x, y))

    return p2

print(Newton(InitialParameters,Xi,Yi))
#print(Newton(InitialParameters,Xi,Yi)[1][1])

#print(ExtremumCondition(InitialParameters,Xi,Yi))
#print(Newton(InitialParameters,Xi,Yi)[0])
#print(Newton(InitialParameters,Xi,Yi)[1])
#print(Newton(InitialParameters,Xi,Yi)[2])
#print(Newton(InitialParameters,Xi,Yi))





