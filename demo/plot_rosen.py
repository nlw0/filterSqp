from pylab import *

def rosenbrock(x, y): return (1 - x)**2 + 100 * (y - x**2)**2

ion()
qq = genfromtxt(sys.argv[1], delimiter='\t')

n=1000

x = np.linspace(-2., 2., n)
y = np.linspace(-1., 3., n)
X, Y = np.meshgrid(x, y)
Z = rosenbrock(X, Y)

figure(figsize=(6.2,6))
contour(X, Y, Z, 0.2 * 2.0 ** mgrid[0:12.0] )
plot(qq[:,1], qq[:,2], '-rs', lw=2)
axis('equal')
grid()
title('Rosenbrock minimization by Newtons method')
