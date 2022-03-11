import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
pi=np.pi
def equispaced(n):
    return np.linspace(-1, 1, n)
def chebyshev(n):
    X = np.zeros(n)
    for i in range(0, n):
        X[i] = np.cos((2*i + 1)/(2*n) * pi)
    return X
def size(X):
    return X.size
def chebyshev_poly(n, x):
    mul=n * np.arccos(x)
    return np.cos(mul)
def monomial(n, x):
    return x**n
def vandermonde(X, m):
    n = size(X)
    V_m = np.zeros((size(X), m))
    V_c = np.zeros((size(X), m))
    for i in range(size(X)):
        for j in range(m):
            V_m[i, j] = monomial(j, X[i])
            V_c[i, j] = chebyshev_poly(j, X[i])
    return V_m, V_c

values={}
values["mono_equi"]=[]
values["cheb_qui"]=[]
values["mono_cheb"]=[]
values["cheb_cheb"]=[]

keys=[]
d_value=[]
for n in range(5, 101, 5):
    X_equi = equispaced(n)
    X_cheb = chebyshev(n)
 
    mono_equi, cheb_equi = vandermonde(X_equi, n)
    mono_cheb, cheb_cheb = vandermonde(X_cheb, n)
    

    values["mono_equi"].append([n, la.cond(mono_equi)])
    values["cheb_qui"].append([n, la.cond(cheb_equi)])
    values["mono_cheb"].append([n, la.cond(mono_cheb)])
    values["cheb_cheb"].append([n, la.cond(cheb_cheb)])
keys=values.keys()
np_mono_equi = np.asarray(values["mono_equi"])
np_cheb_equi = np.asarray(values["cheb_qui"])
np_mono_cheb = np.asarray(values["mono_cheb"])
np_cheb_cheb = np.asarray(values["cheb_cheb"])
d_value=values.values()


plot0_mono_eqi=np_mono_equi[:, 0]
plot0_cheb_eqi=np_cheb_equi[:, 0]
plot0_mono_cheb=np_mono_cheb[:, 0]
plot0_cheb_cheb=np_cheb_cheb[:, 0]

plot1_mono_eqi=np_mono_equi[:, 1]
plot1_cheb_eqi=np_cheb_equi[:, 1]
plot1_mono_cheb=np_mono_cheb[:, 1]
plot1_cheb_cheb=np_cheb_cheb[:, 1]

plt.semilogy(plot0_mono_eqi, plot1_mono_eqi, label='Mono,Equi Nodes',color='blue',marker='.')
plt.semilogy(plot0_cheb_eqi, plot1_cheb_eqi, label='Mono,Cheb Nodes',color='green',marker='.')
plt.semilogy(plot0_mono_cheb, plot1_mono_cheb, label='Cheb_poly,Equi Nodes',color='orange',marker='.')
plt.semilogy(plot0_cheb_cheb, plot1_cheb_cheb, label='Cheb_poly,Cheb Nodes',color='purple',marker='.')
plt.xlabel("Number of Interpolation Nodes",color="#641E16")
plt.ylabel('Condition-number',color="#641E16")
plt.title("Generalized Vandermonde Condition-Number",color="Red")
plt.xticks(visible = True)
plt.legend(loc='best')
plt.savefig("problem_3.png")
plt.show()
