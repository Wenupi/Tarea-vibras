"""Este script calcula un modo a la vez, si se quiere otro
plot hay que cambiar "k" y "n" a mano y volver a correrlo
"""
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import jv
from scipy.optimize import fsolve
from mpmath import besseljzero

"""
Constantes
"""
alpha = np.pi/4  # apertura
radio = 1  # membrana de radio 1
c = 0.75  # velocidad del sonido?

# Tiempos
FPS = 10  # cuadros por segundo
tiempo_por_modo = 10
cuadros = tiempo_por_modo*FPS  # cuadros


def m_angular(k):
    """
    Valor de "m" dado por las condiciones de borde sobre la función
    Theta(theta) (Theta(0)=0 ^ Theta(2pi-alpha)=0). Donde "k" es el
    número entero que sale de la segunda condición de borde.
    """
    return k*np.pi/(2*np.pi-alpha)


def omega_kn(k, n):
    """
    Calcula ω igual al n-1ésimo cero de la función dividido por el
    radio y multiplicado por c (es despejar ω)
    de Bessel de orden m.

    """
    return float(besseljzero(m_angular(k), n))*c/radio


def Zeta(r, theta, t, k, n):
    """
    Otorga la expresión de la deformación de la membrana con el método
    de separación de variables, Zeta=R(r)*Theta(θ)*T(t)
    """
    omega = omega_kn(k, n)
    
    T = np.sin(omega*t)  # T(t)
    R = jv(m_angular(k), omega*r/c)  # R(r)
    Theta = np.sin(m_angular(k) * theta)  # Theta(theta)

    return R * T * Theta

"""
Arrays
"""
r = np.linspace(0, radio, 50)  # distancias
theta = np.linspace(0, 2*np.pi-alpha, 50)  # ángulos
# Para generar el plano
r, theta = np.meshgrid(r, theta)
x = np.cos(theta) * r
y = np.sin(theta) * r

"""
Animación
"""
fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
ax.set_axis_off()


def paso_de_tiempo(i):
    """
    Función que permite calcular Zeta para cada elemento del plano, en
    cada tiempo
    """
    print(f'{i / cuadros:.2%}', end='\r')  # indicar la carga 
    t = i / FPS  # tiempo a evaluar
    ax.cla()
    ax.set_axis_off()
    """
    Aquí se define k y n
    ====================
    """
    k = 5
    n = 1
    """
    ====================
    """
    z = Zeta(r, theta, t, k, n)  # calcula la deformación en cada punto
    vmax = np.max(jv(m_angular(k), np.linspace(0, omega_kn(k, n), 100)))
    # Plotea la superficie
    ax.plot_surface(x, y, z,
                    linewidth=0, cmap='Spectral', vmin=-vmax,
                    vmax=vmax, rcount=100, ccount=100,)

    ax.set_zlim(-1.1, 1.1)
    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(-0.75, 0.75)
    omega = omega_kn(k, n)
    ax.set_title(
        f'Membrana pacman, m = {m_angular(k)}, ω={omega:.2f}',
        size=36, weight='bold', family='Fira Sans',
    )

ani = FuncAnimation(fig, paso_de_tiempo, frames=cuadros, interval=1000/FPS, repeat=False)
ani.save(f'k5-n1.gif', writer='ffmpeg')  # ir cambiando el nombre
