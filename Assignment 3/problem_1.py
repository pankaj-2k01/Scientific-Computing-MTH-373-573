import numpy as np
import matplotlib.pyplot as plt
def newton(func,func_no ,x_0, epsilon):
    x_n = x_0
    while func(x_n,'function',func_no) > epsilon:
        y = func(x_n,'function',func_no)
        y_ = func(x_n,'derivative',func_no)
        x_n = x_n-y/y_
    return x_n

def func(x,choice,function):
    if function==1:
        if choice=='function':
            return x**2 - 1
        elif choice=='derivative':
            return 2*x
    elif function ==2:
        if choice=='function':
            return (x - 1)**4
        elif choice=='derivative':
            return 4*(x - 1)**3
    elif function==3:
        if choice=='function':
            return x - np.cos(x)
        elif choice=='derivative':
            return 1 + np.sin(x)


print(newton(func,1,10**6, 1e-100))
print(newton(func,2,10, 1e-60))
print(newton(func,3,1, 1e-100))


