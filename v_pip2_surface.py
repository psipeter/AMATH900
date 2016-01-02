from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

params={
	'z' : 1.866,
	'F' : 96480,
	'R' : 8314,
	'T' : 300,
	'theta' : 10000,
	'kg' : 0.02,
	'kv0' : 4.512,
	'kp0' : 0.75
}

kv0=params['kv0']
z=params['z']
F=params['F']
R=params['R']
T=params['T']
kp0=params['kp0']
kg=params['kg']
theta=params['theta']

v = np.linspace(-120e-3, 20e-3, 10)
pip2 = np.logspace(-3, 3, 10)

print "\nv", v
print "\npip2", pip2

kv=kv0*np.exp(z*F*(v*1e3)/(R*T))
kp = kp0 * pip2

print "\nkv", kv
print "\nkp", kp

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(v, pip2)

Z = (kg+kv*kg+kp*kg+theta*kv*kp*kg)/(1+kg+kp+kv+kv*kg+kp*kg+kv*kp+theta*kv*kp*kg) / (kg+theta*kp*kg)/(1+kg+kp+theta*kp*kg)
print "z", Z

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
