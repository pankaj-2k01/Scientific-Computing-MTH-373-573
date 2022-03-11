import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

max_iteration= 100000
value=1e-12
ith_steps = []

 
def newton_method(function, Jacobian, x0,max_iteration=500):
    prev = x0
    jth_point = Jacobian(prev)
    ith_point = function(prev)
    store=[]
    y = solve(jth_point, ith_point)
    for i in range(1,max_iteration+1):
        store.append(la.norm(next-previous_main))
        if la.norm(next-prev) > value:
            steps.append(prev)
            store.append(la.norm(next-previous_main))
            y = solve(jth_point, ith_point)
           
    print('Spherical Roots:', steps[-1])
    return steps[-1]

def newton_method_2(function, Jacobian, x0,max_iteration=500):
    previous_main = x0
    jth_point = Jacobian(previous_main)
    ith_point = function(previous_main)
    
    y = solve(jth_point, ith_point)
    next = previous_main-y
 
    for i in range(1,max_iteration+1):
        store.append(la.norm(next-previous_main))
        if la.norm(next-previous_main) > value:
            ith_steps.append(previous_main)
            previous_main = next
            y = solve(jth_point, ith_point)
            store.append(la.norm(next-previous_main))
            next = previous_main-y
    r=ith_steps[-1]
    return r

