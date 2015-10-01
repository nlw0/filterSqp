from pylab import *


ion()

x,y = genfromtxt(sys.argv[1]).T

tt = mgrid[0:1.0:0.001]*2*pi

plot(x, y, '-o')
plot(sin(tt), cos(tt), 'k--')
grid()
axis('equal')
