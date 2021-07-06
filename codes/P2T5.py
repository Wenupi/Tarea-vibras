import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

"""
Parámetros
"""
c = 333  # m/s
rho_ambient = 1.2  # kg/m^3
A = 10  # Pa
L = 15  # m
v_x_in = rho_prime_in = 0


def p_in(x):
    array = np.logical_and(-L/2 <= x, x <= L/2)
    return np.where(array, A, 0)


def f(x, t):
    if t==0:
        return p_in(x)/2
    else:
        array = np.logical_and(-L/2 + c*t <= x, x <= c*t + L/2)
        return np.where(array, A, 0)


def g(x, t):
    if t==0:
        return p_in(x)/2
    else:
        array = np.logical_and(-L/2 - c*t <= x, x <= -c*t + L/2)
        return np.where(array, A, 0)


def p(x, t):
    return f(x, t) + g(x, t)


def rho_prime(x, t):
    return p(x, t)/2 + rho_prime_in - p_in(x)/(c**2)


def rho(x, t):
    return rho_prime(x, t) + rho_prime_in


def v_x(x, t):
    return 1/(rho_ambient*c)*(f(x, t) - g(x, t))*10


x_list = np.linspace(-50 - L/2, 50 + L/2, 1000)
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 1
plt.rcParams['grid.color'] = "#cccccc"
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
tiempo = 3*L/(2*c)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(x_list, p(x_list, tiempo), label='p(x) [Pa]')
ax.plot(x_list, rho_prime(x_list, tiempo), label=r'$\rho^\prime$(x) [$kg/{m}^3$]')
ax.plot(x_list, v_x(x_list, tiempo), label='$v_x$(x) [$m/s$]')
ax.set_xlabel('Posición [m]')
ax.set_title('Campos acústicos para t=3L/(2c)[s]')
ax.legend()
fig.tight_layout()
fig.show()
fig.savefig('plotP2T5.png', dpi=1000)
