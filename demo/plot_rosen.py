from pylab import *

def rosenbrock(x, y): return (1 - x)**2 + 100 * (y - x**2)**2
def himmelblau(x, y): return (x**2 + y - 11)**2 + (x + y**2 - 7)**2
def saddle(x,y): return x * x - y * y
def bowl(x,y): return x * x*2 + y*y

ion()

n=1000

cc = 0.5 * 2.0 ** mgrid[0:10.0:0.5]
x = np.linspace(-2., 2., n)
y = np.linspace(-1., 3., n)
X, Y = np.meshgrid(x, y)
Z = rosenbrock(X, Y)

# cc = mgrid[0.2: 500:5]
# x = np.linspace(-4., 4., n)
# y = np.linspace(-4., 4., n)
# X, Y = np.meshgrid(x, y)
# Z = himmelblau(X, Y)

# cc = mgrid[-4:4:0.05]
# x = np.linspace(-2., 2., n)
# y = np.linspace(-2., 2., n)
# X, Y = np.meshgrid(x, y)
# Z = saddle(X, Y)

figure(figsize=(6.2,6))
contour(X, Y, Z, cc )

for filename in sys.argv[1:]:
    qq = genfromtxt(filename, delimiter='\t')
    plot(qq[:,4], qq[:,5], '-ro', lw=2)
    plot(qq[0:1,4], qq[0:1,5], 'bs', lw=2)
    plot(qq[-1:,4], qq[-1:,5], 'ys', lw=2)

axis('equal')
grid()
xlim(x[0], x[-1])
ylim(y[0], y[-1])

method = filename[:-4]

title('Rosenbrock minimization - {}'.format(method))
# title('Himmelblau minimization - fixed trust region')
# title('Minimization - {}'.format(method))

savefig(method + ".png")
