# Problem 2: Gaussian elimination without and with partial
# pivoting.

import numpy as np
from scipy.linalg import hilbert

# Gaussian Elimination without partial pivoting for factorizing a
# linear system A x = b into A = L * U and b = L^{-1}x * b
#
# Inputs: Matrix A, Vector b
#
# Outputs: Solution x
#
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


# Gaussian Elimination with partial pivoting for factorizing a
# linear system A x = b into P * A = L * U and b = L^{-1}x * P * b
#
# Inputs: Matrix A, Vector b
#
# Outputs: Solution x
#
# Note: The permutation matrix is not tracked
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

n = 10
matrix_choice = 1   # 1, 2 or 3. Default: random

def matrix_choice(choice ,n):
    A=np.array([])
    if choice == 2:
        A=hilbert(n)
    elif choice == 3:
        A=np.ones((n,n))
        for i in range(n):
            for j in range(i):
                A[i][j]*=-1

    else:
        # Random matrix example
        np.random.seed(0)
        A = np.random.rand(n, n)

    return A

L=[i for i in range(10,50,10)]
# for choice in range(1,4):
print("----------------------------------------------------------------------------")
for n in L:
    print("\t\t\t\t\t\t\tWhen N:",n)
    print("----------------------------------------------------------------------------")
    for choice in range(1,4):
        A=matrix_choice(choice,n)

        x_star = np.ones(n)
        b = np.dot(A, x_star)
            
        # Computations as sought
        
        print("choice :",choice)
        
        print ("condition number: %g" % np.linalg.cond(A))
        

        # Solve without pivoting
        x_npp = GE(A.copy(), b.copy())
        error_npp = np.linalg.norm(x_star - x_npp)
        residual_npp = np.linalg.norm(np.dot(A, x_npp) - b)
        print ("No partial pivoting:","Error = %1.6g  Residual = %1.6g" %(error_npp, residual_npp))
       
        #print ("Error = %1.6g  Residual = %1.6g" %(error_npp, residual_npp))

        # Solve with partial pivoting
        x_pp = GE_pp(A.copy(), b.copy())
        error_pp = np.linalg.norm(x_star - x_pp)
        residual_pp = np.linalg.norm(np.dot(A, x_pp) - b)
        print ("With partial pivoting:","Error = %1.6g  Residual = %1.6g" %(error_pp, residual_pp))
        
        #print ("Error = %1.6g  Residual = %1.6g" %(error_pp, residual_pp))

        # Solve using numpy.linalg's solve
        x_lin=np.linalg.solve(A.copy(),b.copy())
        error_lin = np.linalg.norm(x_star - x_lin)
        residual_lin = np.linalg.norm(np.dot(A, x_lin) - b)
        print ("Using linalg:","Error = %1.6g  Residual = %1.6g" %(error_lin, residual_lin))
        print("----------------------------------------------------------------------------")

        #print ("Error = %1.6g  Residual = %1.6g" %(error_npp, residual_npp))


