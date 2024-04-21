import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from sympy import symbols, Eq, solve

plt.close('all')

# chosen coefficients
a = 1.1        # prey population growth rate
b = 0.4        # prey population death rate
c = 0.1        # predator population growth rate
d = 0.4        # predator population death rate

# lotka_volterra function
def lotka_volterra_ivp(t, Y):
    dxdt = (a - b * Y[1]) * Y[0]
    dydt = (c*Y[0] - d) * Y[1]
    return np.array([dxdt, dydt])

def lotka_volterra_odeint(Y, t):
    x, y = Y
    dxdt = (a - b*y) * x
    dydt = (c*x - d) * y
    return [dxdt, dydt]

# initial conditions
initial_sizes = np.array([10, 10])
t0 = 0
tk = 100
Nt = 1000
t = np.linspace(t0, tk, Nt)

# solving the ODE
sol_ivp = solve_ivp(lotka_volterra_ivp, (t0, tk), initial_sizes, t_eval=t, atol=1e-8, rtol=1e-8)
sol_odeint = odeint(lotka_volterra_odeint, initial_sizes, t)

# plotting the solution
plt.figure(figsize=(10, 25))
# plotting the populations vs time 
plt.plot(sol_ivp.t, sol_ivp.y[0], label='prey population')
plt.plot(sol_ivp.t, sol_ivp.y[1], label='predator population')
plt.xlabel('time')
plt.ylabel('population')
plt.title('Lotka-Volterra model using solve_ivp')
plt.legend()
plt.grid(True)

plt.figure(figsize=(10, 25))
plt.plot(t, sol_odeint[:, 0], label='prey population')
plt.plot(t, sol_odeint[:, 1], label='predator population')
plt.xlabel('time')
plt.ylabel('population')
plt.title('Lotka-Volterra model using odeint')
plt.legend()
plt.grid(True)
plt.show()

# plotting the phase plot (populations against each other - predator(prey))
plt.figure(figsize=(10, 10))
plt.plot(sol_ivp.y[0], sol_ivp.y[1])
plt.xlabel('prey population')
plt.ylabel('predator population')
plt.title('phase plot using solve_ivp')
plt.grid(True)

plt.figure(figsize=(10, 10))
plt.plot(sol_odeint[:, 0], sol_odeint[:, 1])
plt.xlabel('prey population')
plt.ylabel('predator population')
plt.title('phase plot using odeint')
plt.grid(True)
plt.show()

# checking the equilibrium points

# define symbols
x, y = symbols('x y')

# lotka-volterra equations equaled to zero
equation1 = Eq((a - b*y)*x, 0)
equation2 = Eq((c*x - d)*y, 0)

# solving for equilibrium points with the use of 'solve' command 
sol = solve((equation1, equation2), (x, y))
equilibrium_points = (sol[-1])
print("equilibrium points:", equilibrium_points)