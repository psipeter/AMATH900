import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def randrange(n, vmin, vmax):
    return (vmax - vmin)*np.random.rand(n) + vmin

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
n = 100
# for c, m, zl, zh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
#     xs = randrange(n, 23, 32)
#     ys = randrange(n, 0, 100)
#     zs = randrange(n, zl, zh)
#     ax.scatter(xs, ys, zs, c=c, marker=m)

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

v = np.linspace(-120e-3, 80e-3, 40)
pip2 = np.logspace(-6, 3, 40, base=10.0)
# print "\nv", v
# print "\npip2", pip2
kv = kv0*np.exp(z*F*(v*1e3)/(R*T))
kp = kp0 * pip2

for i in range(len(kv)):
	for j in range(len(kp)):
		xs = v[i]
		ys = np.log10(pip2[j])
		zs = (((kg+kv[i]*kg+kp[j]*kg+theta*kv[i]*kp[j]*kg) / 
			(1+kg+kp[j]+kv[i]+kv[i]*kg+kp[j]*kg+kv[i]*kp[j]+theta*kv[i]*kp[j]*kg))
			/ ((kg+theta*kp[j]*kg)/(1+kg+kp[j]+theta*kp[j]*kg)))
		# print xs
		# print ys
		# print zs
		ax.scatter(xs, ys, zs)

ax.set_xlabel('voltage')
ax.set_ylabel('log_10 kcnq_pip2')
ax.set_zlabel('kcnq channel open')

plt.show()
