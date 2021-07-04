import numpy as np
import matplotlib.pyplot as plt

c = 333  # m/s
rho_ambient = 1.2  # kg/m^3
A = 10
L = 10
v_x_in = rho_prime_in = 0


def p_in(x):
    array = np.logical_and(-L/2 <= x, x <= L/2)
    return np.where(array, A, 0)
    #if (-L/2 <= x) and (x <= L/2):
    #    return A
    #else:
    #    return 0


def f(x, t):
    #return (p_in(x)*(x - c*t) + rho(x, t)*c*v_x_in*(x - c*t))/2
    return (p_in(x) + 0)/2


def g(x, t):
    #return (p_in(x)*(x + c*t) + rho(x, t)*c*v_x_in*(x + c*t))/2
    return (p_in(x) + 0)/2


def p(x, t):
    return f(x, t) + g(x, t)


def rho_prime(x, t):
    return p(x, t)/2 + rho_prime_in - p_in(x)/(c**2)


def rho(x, t):
    return rho_prime(x, t) + rho_prime_in


def v_x(x, t):
    return 1/(rho_ambient*c)*(f(x, t) - g(x, t))


x_list = np.linspace(-10 - L/2, 10 + L/2, 1000)
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 1
plt.rcParams['grid.color'] = "#cccccc"
plt.rcParams.update({"font.size": 9, "font.family": "serif"})

#p_array = []
#rho_prime_array = []
#v_x_array = []
#for x in x_list:
#    p_array.append(p(x, 3*L/(2*c))) 
#    rho_prime_array.append(rho_prime(x, 3*L/(2*c)))
#    v_x_array.append(v_x(x, 3*L/(2*c)))

fig, ax = plt.subplots(figsize=(5, 3.5))
ax.plot(x_list, p(x_list, 3*L/(2*c)), label='Presión')
#ax.plot(x_list, p_array, label='Presión')
ax.plot(x_list, rho_prime(x_list, 3*L/(2*c)), label='Densidad prima')
#ax.plot(x_list, rho_prime_array, label='Densidad prima')
ax.plot(x_list, v_x(x_list, 3*L/(2*c)), label='Velocidad')
#ax.plot(x_list, v_x_array, label='Velocidad')
ax.set_xlabel('Posición [m]')
ax.legend()
fig.show()