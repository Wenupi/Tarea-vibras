import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import jv
from scipy.optimize import fsolve

"""
Constantes
"""
alpha = np.pi/4  # apertura
radio = 1  # membrana de radio 1
c = 0.75  # velocidad del sonido?
# Este "m" debería ir multiplicado por un k
# pertenenciente a los enteros. Lo definí como k=1
m = np.pi/(2*np.pi-alpha)
"""
Se supone que debería haber una función que calcule
los ceros de J_n, pero solo acepta n entero, así
que en el resto del script se ocupó un único cero
sacado de tabla
"""
# Ceros de la función de Bessel de orden "m"
#ceros_J_m = [jn_zeros(m, 10) for m in range(10)]
#ceros_J_solo_m = jn_zeros(m, 10)
primer0 = 3.242444822093764691446
segundo0 = 0
tercer0 = 9.532984858265578140957

# Tiempos
FPS = 10  # cuadros por segundo
tiempo_por_modo = 10
cuadros = tiempo_por_modo*FPS  # cuadros

# ESTA FUNCIÓN NO FUNCIONA HASTA ENCONTRAR UNA FUNCIÓN QUE CALCULE CEROS DE J_n
def omega_mn(m, n):
    """
    Calcula omega'=omega/c, donde radio*omega' es igual al n-1ésimo cero de la función
    de Bessel de orden m.
    *Nota mental:
    - Si escojo solo un cero de J_m voy a tener solo un modo de oscilación,
    cambiar el n permite calcular las otras frecuencias de oscilación. (parece que no es cierto)
    - Omega depende tanto de m (orden de la función de Bessel) como del cero escogido para esa
    esa función
    - Este "n" no es el mismo "n" que tengo en la tarea
    """
    return ceros_J_solo_m[n - 1]*c/radio


def Zeta(r, theta, t):
    """
    Otorga la expresión de la deformación de la membrana con el método
    de separación de variables, Zeta=R(r)*Omega(omega)*T(t)
    """
    # Omega calculado con tablas para el "m" dado 
    omega_prima = tercer0*c/radio  # puede ser primer0, segundo0 o tercer0
    
    T = np.sin(omega_prima*t)  # T(t)
    R = jv(m, omega_prima*r/c)  # R(r)
    Theta = np.sin(m * theta)  # Theta(theta)

    return R * T * Theta

"""
Arreglos
"""
r = np.linspace(0, radio, 50)  # distancias
theta = np.linspace(0, 2*np.pi-alpha, 50)  # ángulos
# Para generar el plano (proyección?)
r, theta = np.meshgrid(r, theta)
x = np.cos(theta) * r
y = np.sin(theta) * r

"""
Animación
"""
fig = plt.figure(figsize=(19.2, 10.8), dpi=100)
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
ax.set_axis_off()


# La función a continuación también debería permitir ir cambiando
# de modo, pero no hasta tener la función de ceros para J_n

def paso_de_tiempo(i):
    """
    Función que permite calcular Zeta para cada elemento del plano, en
    cada tiempo
    """
    print(f'{i / cuadros:.2%}', end='\r')  # indicar la carga 
    t = i / FPS  # tiempo a evaluar
    ax.cla()
    ax.set_axis_off()
    z = Zeta(r, theta, t)  # calcula la deformación en cada punto
    vmax = np.max(jv(m, np.linspace(0, tercer0, 100)))
    # Plotea la superficie
    ax.plot_surface(x, y, z,
                    linewidth=0, cmap='Spectral', vmin=-vmax,
                    vmax=vmax, rcount=100, ccount=100,)

    ax.set_zlim(-1.1, 1.1)
    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(-0.75, 0.75)
    omega = tercer0
    ax.set_title(
        f'Membrana pacman, m = {m}, ω={omega:.2f}',
        size=36, weight='bold', family='Fira Sans',
    )

ani = FuncAnimation(fig, paso_de_tiempo, frames=cuadros, interval=1000/FPS, repeat=False)
ani.save(f'membrana_pacman3.gif', writer='ffmpeg')
