from pylab import *

def thefun(alpha, lam, nu):
    return norm(alpha / (lam + nu))

ion()

input = open(sys.argv[1]).readlines()

rho = float(input[0])
n = int(input[1])
alpha = array([float(x) for x in input[2].split(' ')])
lam = array([float(x) for x in input[3].split(' ')])

nu = mgrid[max(-min(lam), 0)+0.01:100:0.01]

yy = array([thefun(alpha, lam, n) for n in nu])

semilogy(nu, yy, lw=2)
semilogy(nu[[0, -1]], [rho, rho], 'r--', lw=2)
grid()
