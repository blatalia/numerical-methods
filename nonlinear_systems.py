"""
authors:
task: 7th task - solving nonlinear systems of equations
equations:
---first Kirchhoff law---
I1 = I6 + I5
I3 = I1 + I2
I4 = I2 + I5
---second Kirchhoff law---
E4 = I3*R6 + I4*R5 + I2*r4
E1 = I3*R6 + I6*R7 + I1*r1
E1 - E4 = I5*R3 + I5*R2 + I1*r1 - I2*r4
"""
import numpy as np
from scipy.optimize import fsolve

def kirchhoff_equations(X):
    [I1, I2, I3, I4, I5, I6] = X

    eq1 = I1 - I6 - I5 
    eq2 = I3 - I1 - I2 
    eq3 = I4 - I2 - I5

    eq4 = E4 - I3 * R6_0 - I4 * R5_0 - I2 * r4_0 
    eq5 = E1 - I3 * R6_0 - I6 * R7_0 - I1 * r1_0 
    eq6 = E1 - E4 - I5 * R3_0 - I5 * R2_0 - I1 * r1_0 + I2 * r4_0 

    return np.array([eq1, eq2, eq3, eq4, eq5, eq6])   

# values of resistances [Ohm]: 
r1_0 = 4 
r4_0 = 5 
R2_0 = 14 
R3_0 = 7 
R5_0 = 17 
R6_0 = 11 
R7_0 = 20
# values of voltage [V]: 
E4 = 30
E1 = 22 

starting = np.array([0, 0, 0, 0, 0, 0])

sol = fsolve(kirchhoff_equations, starting, xtol=1e-9, maxfev=1000)

I1_0 = sol[0]
I2_0 = sol[1]
I3_0 = sol[2]
I4_0 = sol[3]
I5_0 = sol[4]
I6_0 = sol[5]

for _ in range(len(sol)):
    print(f"the value of current I{_+1} is {sol[_]}")
print("\n")

# updated values of resistances [Ohm]:
kp = 0.1
kn = -0.1

r1_1 = r1_0 + I1_0 * kn
r4_1 = r4_0 + I4_0 * kn
R2_1 = R2_0 + I2_0 * kp
R3_1 = R3_0 + I3_0 * kn
R5_1 = R5_0 + I5_0 * kn
R6_1 = R6_0 + I6_0 * kp
R7_1 = R7_0 + I6_0 * kp

def kirchhoff_updated(Y):
    [I1_1, I2_1, I3_1, I4_1, I5_1, I6_1] = Y

    eq1 = I1_1 - I6_1 - I5_1
    eq2 = I3_1 - I1_1 - I2_1
    eq3 = I4_1 - I2_1 - I5_1

    eq4 = E4 - I3_1 * R6_1 - I4_1 * R5_1 - I2_1 * r4_1
    eq5 = E1 - I3_1 * R6_1 - I6_1 * R7_1 - I1_1 * r1_1
    eq6 = E1 - E4 - I5_1 * R3_1 - I5_1 * R2_1 - I1_1 * r1_1 + I2_1 * r4_1

    return np.array([eq1, eq2, eq3, eq4, eq5, eq6])

sol_updated = fsolve(kirchhoff_updated, starting, xtol=1e-9, maxfev=1000)

for e in range(len(sol_updated)):
    print(f"the new value of current I{e+1} is {sol_updated[e]}")
print("\n")

def comparison(first_sol,second_sol):
    comparison = []
    for i in range(len(first_sol)):
        comparison.append(abs(first_sol[i] - second_sol[i]))
    return comparison
comp = comparison(sol, sol_updated)

for i in range(len(comp)):
    print(f"for I{i+1} the difference between solutions is {comp[i]}")
print("\n")

check_part_a = kirchhoff_equations(sol)
print(check_part_a)
check_part_b = kirchhoff_updated(sol_updated)
print(check_part_b)


