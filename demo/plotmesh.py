from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


ion()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x,y,z = genfromtxt(sys.argv[1]).T

x = x.reshape(11, -1)
y = y.reshape(11, -1)
z = z.reshape(11, -1)


ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=1, antialiased=False)

ax = fig.gca(projection='3d')
ax.set_aspect('equal')
MAX = 3
for direction in (-1, 1):
    for point in np.diag(direction * MAX * np.array([1,1,1])):
        ax.plot([point[0]], [point[1]], [point[2]], 'w')
