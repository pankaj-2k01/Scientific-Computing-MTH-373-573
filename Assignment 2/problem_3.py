import numpy as np
from scipy.linalg import hilbert

def GE(A, b):
    n=len(A)
    x=np.ones(n)
    for k in range(n-1):
        for i in range(k+1,n):
            temp=A[i][k]/A[k][k]
            A[i][k]=temp
            for j in range(k+1,n):
                A[i][j]=A[i][j]-(temp*A[k][j])
            b[i]=b[i]-(temp*b[k])

    x[n-1]=b[n-1]/A[n-1][n-1]
    for i in range(n-2,-1,-1):
        sum=b[i]
        for j in range(i+1,n):
            sum=sum-(A[i][j]*x[j])
        x[i]=sum/A[i][i]
    
    
    return x
def GE_pp(A, b):
    n=len(A)
    
    
    x=np.ones(n)
    l=np.ones(n,dtype=int)
    s=np.ones(n)
    for i in range(n):
        l[i]=i
        smax=0
        for j in range(n):
            smax=max(smax,abs(A[i][j]))
        s[i]=smax

    for k in range(n-1):
        rmax=0
        index=0
        for i in range(k,n):
            r=abs((A[l[i]][k])/s[l[i]])
            if r>rmax:
                rmax=r
                index=i
        ltemp=l[k]
        l[k]=l[index]
        l[index]=ltemp
        for i in range(k+1,n):
            amult=((A[l[i]][k])/(A[l[k]][k]))
            A[l[i]][k]=amult
            for j in range(k+1,n):
                A[l[i]][j]=A[l[i]][j]-(amult*A[l[k]][j])
    for k in range(n-1):
        for i in range(k+1,n):
            b[l[i]]=b[l[i]]-(A[l[i]][k]*b[l[k]])
    x[n-1]=b[l[n-1]]/A[l[n-1]][n-1]
    for i in range(n-2,-1,-1):
        sum=b[l[i]]
        for j in range(i+1,n):
            sum=sum-(A[l[i]][j]*x[j])
        x[i]=sum/A[l[i]][i]

    return x
#-----------------Q3a
A=np.array([[1,1,2e9],[2,-1,1e9],[1,2,0]])
B=np.array([1.0,1.0,1.0])
print(GE(A,B))

#------------------Q3b
A_1=np.array([[1,1,2e9],[2,-1,1e9],[1,2,0]])
B_1=np.array([1.0,1.0,1.0])
for i in range(len(A_1)):
    maxValue=abs(max(A_1[i],key=abs))
    B_1[i]=B_1[i]/maxValue
    A_1[i]=A_1[i]/maxValue
print(GE_pp(A_1,B_1))