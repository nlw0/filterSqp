from pylab import *

def rosenbrock(x, y): return (1 - x)**2 + 100 * (y - x**2)**2
def himmelblau(x, y): return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

ion()
filename = sys.argv[1]
qq = genfromtxt(filename, delimiter='\t')

n=1000

# x = np.linspace(-2., 2., n)
# y = np.linspace(-1., 3., n)
# X, Y = np.meshgrid(x, y)
# Z = rosenbrock(X, Y)

x = np.linspace(-4., 4., n)
y = np.linspace(-4., 4., n)
X, Y = np.meshgrid(x, y)
Z = himmelblau(X, Y)

figure(figsize=(6.2,6))
#cc = 0.5 * 2.0 ** mgrid[0:10.0:0.1]
cc = mgrid[0.2: 500:5]
#cc = [0.2, 1.2, 10.2, 20.2, 120, 250, 500]
contour(X, Y, Z, cc )
plot(qq[:,1], qq[:,2], '-rs', lw=2)
axis('equal')
grid()
xlim(x[0], x[-1])
ylim(y[0], y[-1])

method = filename[:-4]

# title('Rosenbrock minimization - {}'.format(method))
title('Himmelblau minimization - {}'.format(method))

savefig(method + ".png")
