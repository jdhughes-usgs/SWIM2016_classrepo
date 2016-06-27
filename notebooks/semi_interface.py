from __future__ import division
from __future__ import print_function
from pylab import *
from scipy.special import ellipkinc as F
from scipy.special import ellipeinc as E
from scipy.optimize import fsolve, brentq


def fint(phi, a):
    kappasq = (2 + sqrt(3)) / 4  # Eq. 43
    kappa = sqrt(kappasq)
    # theta = arccos((sqrt(3) - 1 - phi / a) / (sqrt(3) + 1 + phi / a))  # Eq. 43
    theta = arccos(
        -1 + 2 * sqrt(3) / (sqrt(3) + 1 + phi / a))  # Eq. 43, reworked
    # Eq. 42
    rv = (3 ** (-0.25) - 3 ** 0.25) * F(theta, kappasq) + \
         2 * 3 ** 0.25 * E(theta, kappasq) - \
         2 * 3 ** 0.25 * sin(theta) * sqrt(1 - kappasq * sin(theta) ** 2) / (
         1 + cos(theta))
    return rv


class SemiCoast:
    def __init__(self, k, H, c, grad, rhof, rhos, L, ztop=0, sealevel=0):
        self.k = k
        self.H = H
        self.c = c
        self.grad = grad
        self.rhof = rhof
        self.rhos = rhos
        self.L = L
        self.ztop = ztop
        self.hs = sealevel + (sealevel - ztop) * (self.rhos - self.rhof) / self.rhof
        self.initialize()

    def initialize(self):
        # Variables
        self.lab = sqrt(self.k * self.H * self.c)
        self.nu = (self.rhos - self.rhof) / self.rhof
        self.mu = self.grad * self.lab / self.H / self.nu

        # Case I or Case II
        if self.mu < sqrt(2 / 3):
            self.Loverlab = (18 * self.mu) ** (1 / 3)  # Eq. 25
            if self.Loverlab * self.lab <= self.L:
                self.case = 1
                self.doverlab = (1 - (3 * self.mu ** 2 / 2) ** (2 / 3)) / (
                2 * self.mu)  # Eq.28
                self.phi0 = (3 * self.mu ** 2 / 2) ** (1 / 3)
            else:
                self.case = None
        else:
            self.Loverlab = sqrt(6)  # Eq. 29
            self.doverlab = log((self.mu + sqrt(self.mu ** 2 + 1 / 3)) / (
            1 + sqrt(2 / 3)))  # Eq.38
            if self.Loverlab * self.lab + self.doverlab * self.lab <= self.L:
                self.case = 2
                self.phicoast = (1 - sqrt(2 / 3)) / 2 * exp(-self.doverlab) + \
                                (1 + sqrt(2 / 3)) / 2 * exp(self.doverlab)
            else:
                self.case = None

        if self.case is None:
            self.Lsoverlab = self.L / self.lab
            atr = fsolve(self.find_atr, 0.1)  # Eq. 46
            mutr = sqrt(2 * (1 + atr ** 3) / 3)  # Eq. 47
            if self.mu < mutr:
                self.case = 3
                amax = (3 * self.mu ** 2 / 2) ** (1 / 3)
                rv = fsolve(self.find_a_case3, amax / 2, full_output=1)
                if rv[2] == 1:
                    self.a = rv[0]
                    self.phi0 = (3 * self.mu ** 2 / 2 - self.a ** 3) ** (
                    1 / 3)  # Eq. 48
                    #self.doverlab = (1 - (3 * self.mu ** 2 / 2) ** (2 / 3)) / (
                    #2 * self.mu)  # Eq.28
                    self.doverlab = (1 - self.phi0 ** 2) / (
                    2 * self.mu)  # Eq.27 solved for u with phi=1 and phi0
                else:
                    print('Error: a was not found for case 3, mu was:' + str(self.mu))
            else:
                self.case = 4
                rv = brentq(self.find_doverlab_case4, 0, self.Lsoverlab,
                            full_output=1)
                if rv[1].converged == 1:
                    self.doverlab = rv[0]
                    gamma0 = -tanh(self.doverlab) + self.mu / cosh(
                        self.doverlab)  # Eq.53
                    self.a = (3 * gamma0 ** 2 / 2 - 1) ** (1 / 3)  # Eq.54
                    self.phicoast = 1.0 / cosh(self.doverlab) - self.mu * \
                                    sinh(-self.doverlab) / cosh(self.doverlab)
                else:
                    print('Error: doverlab was not found for case 4')

    def toe(self):
        if self.case == 1:
            rv = -self.doverlab * self.lab
        elif self.case == 2:
            rv = self.doverlab * self.lab
        elif self.case == 3:
            rv = -self.doverlab * self.lab
        elif self.case == 4:
            rv = self.doverlab * self.lab
        return rv

    def tip(self):
        if self.case == 1:
            rv = self.Loverlab * self.lab
        elif self.case == 2:
            rv = (self.doverlab + self.Loverlab) * self.lab
        elif self.case == 3:
            rv = self.Lsoverlab * self.lab
        elif self.case == 4:
            rv = self.Lsoverlab * self.lab
        return rv

    def interface(self, N=100):
        # returns (x, z) where z is interface elevation for N points
        if self.case == 1:
            phi = nan * ones(N)
            u = np.linspace(-self.doverlab, self.Loverlab, N)
            u1 = (-self.doverlab <= u) & (u <= 0)
            phi[u1] = sqrt(-2 * self.mu * u[u1] + self.phi0 ** 2)  # Eq. 27
            u2 = (0 < u) & (u <= self.Loverlab)
            phi[u2] = (u[u2] - self.Loverlab) ** 2 / 6
            x = u * self.lab
        elif self.case == 2:
            u = np.linspace(0, self.Loverlab, N)
            phi = (u - self.Loverlab) ** 2 / 6
            x = (self.doverlab + u) * self.lab
        elif self.case == 3:
            phi = linspace(1, 0, N)
            u = nan * ones(N)
            phi1 = phi >= self.phi0
            u[phi1] = (self.phi0 ** 2 - phi[phi1] ** 2) / (2 * self.mu)
            phi2 = phi < self.phi0
            u[phi2] = sqrt(3 * self.a / 2) * (
            fint(phi[phi2], self.a) - fint(0, self.a)) + self.Lsoverlab
            x = u * self.lab
        elif self.case == 4:
            phi = linspace(1, 0, N)
            u = sqrt(3 * self.a / 2) * (fint(phi, self.a) - fint(0, self.a)) + (
            self.Lsoverlab - self.doverlab)
            x = (self.doverlab + u) * self.lab
        h = self.nu * self.H * phi + self.hs
        eta = (h - self.hs) / self.nu
        zeta = self.ztop - eta
        return x, zeta

    def head(self, x):
        # Only implemented for Case 1
        if self.case == 1:
            u = x / self.lab
            if u <= -self.doverlab:
                phi = self.grad * (-self.doverlab - u) + 1
            elif (u > -self.doverlab) & (u <= 0):
                phi = sqrt(-2 * self.mu * u + self.phi0 ** 2)
            elif (u > 0) & (u <= self.Loverlab):
                phi = (u - self.Loverlab) ** 2 / 6
            else:
                phi = nan
        elif self.case == 2:
            pass
        return self.nu * self.H * phi + self.hs

    def onshorex(self, h):
        # returns onshore x location of given head
        phi = (h - self.hs) / (self.nu * self.H)
        if self.case == 1:
            if phi >= 1:
                u = (1 - phi) * (self.nu * self.H) / (
                self.grad * self.lab) - self.doverlab
            elif (phi < 1) & (phi > self.phi0):
                u = (self.phi0 ** 2 - phi ** 2) / (2 * self.mu)
            else:
                u = 0
            x = u * self.lab
        elif self.case == 2:
            if phi > self.phicoast:
                u = (self.phicoast - phi) * (self.nu * self.H) / (
                    self.grad * self.lab) - self.doverlab
            else:
                u = -doverlab
            x = (self.doverlab + u) * self.lab
        elif self.case == 3:
            if phi >= 1:
                u = (1 - phi) * (self.nu * self.H) / (
                self.grad * self.lab) - self.doverlab
            elif (phi < 1) & (phi > self.phi0):
                u = (self.phi0 ** 2 - phi ** 2) / (2 * self.mu)
            else:
                u = 0
            x = u * self.lab
        elif self.case == 4:
            if phi > self.phicoast:
                u = (self.phicoast - phi) * (self.nu * self.H) / (
                    self.grad * self.lab) - self.doverlab
            else:
                u = -self.doverlab
            x = (self.doverlab + u) * self.lab            
        return x

    def find_atr(self, a):
        # Eq. 46. atr is the zero of this function
        rv = sqrt(3 * a / 2) * (fint(1, a) - fint(0, a)) + self.Lsoverlab
        return rv

    def find_a_case3(self, a):
        phi0 = (3 * self.mu ** 2 / 2 - a ** 3) ** (1 / 3)
        rv = sqrt(3 * a / 2) * (fint(phi0, a) - fint(0, a)) + self.Lsoverlab
        return rv

    def find_doverlab_case4(self, doverlab):
        # Function that has to become zero for case IV
        gamma0 = -tanh(doverlab) + self.mu / cosh(doverlab)  # Eq.53
        term = 3 * gamma0 ** 2 / 2 - 1  # Eq.54 (argument only)
        # print('doverlab, term', doverlab, term)
        if term < 0:
            # print('term < 0 in find_doverlab_case4')
            rv = -1
        else:
            a = term ** (1 / 3)  # Eq.54
            rv = sqrt(3 * a / 2) * (
            fint(1, a) - fint(0, a)) + self.Lsoverlab - doverlab
        return rv
    
    def plot(self, xmin=None, xmax=None, newfig=True):
        if newfig:
            if xmin is None:
                xmin = self.toe()
            if xmax is None:
                xmax = self.tip()
            figure(figsize=(8, 4))
            plot([xmin, xmax], [self.ztop - self.H, self.ztop - self.H], linewidth=5, color='k')
            plot([min(xmin, 0), min(xmax, 0)], [self.ztop, self.ztop], linewidth=5, color='k')
            plot([max(xmin, 0), min(xmax, self.L)], [self.ztop, self.ztop], linewidth=5, color=[.8, .8, .8])
            xlim(xmin, xmax)
        x, z = self.interface()
        plot(x, z, zorder=100)
        show()
    
class SemiCoastHead(SemiCoast):
    def __init__(self, k, H, c, h, x, rhof, rhos, L, ztop=0, sealevel=0):
        assert x <=0, "Input error: x must be less than zero"
        self.Linput = L
        # First find position without finite L, then use the specified L
        # This is a bit clunky, but may work ok
        SemiCoast.__init__(self, k, H, c, 0.001, rhof, rhos, inf, ztop, sealevel)
        assert h > self.hs, "Input error: inland head smaller than equivalent freshwater head at top of aquifer"
        self.givenx = x
        self.givenh = h
        # find gradient using log transform of grad to avoid negative values
        start = np.log((self.givenh - self.hs) / abs(x))
        result = fsolve(self.findgrad, start, full_output=1)
        if result[2] == 1:
            self.grad = np.exp(result[0])
            self.initialize()
        else:
            print('Error: gradient could not be found with fsolve when L=inf')
        #
        if self.tip() > self.Linput:
            self.L = self.Linput
            self.initialize()
            result = fsolve(self.findgrad, np.log(self.grad), full_output=1)
            if result[2] == 1:
                self.grad = np.exp(result[0])
                self.initialize()
            else:
                print('Error: gradient could not be found with fsolve when L=inf')
            
    def findgrad(self, g):
        # g is log transformed
        self.grad = np.exp(g)
        self.initialize()
        xh = self.onshorex(self.givenh)
        return xh - self.givenx


## Input
## Case 1
sc1 = SemiCoast(k=10, H=10, c=100, grad=0.0005, rhof=1000, rhos=1025, L=1000)
## Case 2
#sc2 = SemiCoast(k=10, H=10, c=100, grad=0.00375, rhof=1000, rhos=1025, L=1000)
## Case 3
#sc3 = SemiCoast(k=10, H=10, c=100, grad=0.0005, rhof=1000, rhos=1025, L=80)
## Case 4
#sc4 = SemiCoast(k=10, H=10, c=100, grad=0.00375, rhof=1000, rhos=1025, L=150)
#
#def findgrad(g, k=10, H=10, c=100, rhof=1000, rhos=1025, x=-1000, h=1, L=inf):
#    assert x <=0, "Input error: x must be less than zero"
#    sc = SemiCoast(k=k, H=H, c=c, grad=g, rhof=rhof, rhos=rhos, L=L)
#    xh = sc.onshorex(h)
#    return xh - x
#
#g = fsolve(findgrad, 0.001, args=(10, 10, 100, 1000, 1025, -1000, 1, inf))
#
sch = SemiCoastHead(k=10, H=10, c=100, x=-1000, h=1, rhof=1000, rhos=1025, L=1000, ztop=-10, sealevel=0)

    
