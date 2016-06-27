from __future__ import division
from __future__ import print_function
from pylab import *
from scipy.optimize import brentq

class WellNearCoast:
    def __init__(self, k=10, D=20, rhof=1000, rhos=1025, Q0=0.2, d=2000, rw=0.3, Q=800):
        self.k = k
        self.D = D
        self.alpha = rhof / (rhos - rhof)
        self.Q0 = Q0
        self.d = d
        self.rw = rw
        self.Q = Q
        self.phitoe = 0.5 * self.k * (self.alpha + 1) / (self.alpha ** 2) * self.D ** 2
        self.htoe = sqrt(2 * self.phitoe / (self.k * (self.alpha + 1)))
        self.C = -0.5 * k * (self.alpha + 1) / self.alpha * self.D ** 2
        self.xmax = brentq(self.Qxline, 0, self.d - 0.1)
        self.phimax = self.Q / (2 * pi) * log(sqrt((self.xmax - self.d) ** 2) /
                                              sqrt((self.xmax + self.d) ** 2)) + self.Q0 * self.xmax
        self.hmax = sqrt(2 * self.phimax / (self.k * (self.alpha + 1)))
        self.zetamax = -self.alpha * self.hmax
        if self.zetamax > -self.D:
            self.pumpfresh = False
        else:
            self.pumpfresh = True

    def Qxline(self, x):
        # Qx along y=0
        Qx = -self.Q / (2 * pi) * (1 / (x - self.d) - 1 / (x + self.d)) - self.Q0
        return Qx
        # Determine if interface has reached well

    def head_interface_well(self, x, y):
        phi = self.Q / (2 * pi) * log(sqrt((x - self.d) ** 2 + y ** 2) / sqrt((x + self.d) ** 2 + y ** 2)) + self.Q0 * x
        h =  nan * ones_like(x)
        if self.pumpfresh:
            iwel = count_nonzero(x[0] < self.xmax)
            h1 = h[:, :iwel]
            phi1 = phi[:, :iwel]
            h2 = h[:, iwel:]
            phi2 = phi[:, iwel:]
            h1[phi1 <= self.phitoe] = sqrt(2 * phi1[phi1 <= self.phitoe] / (self.k * (self.alpha + 1)))
            h1[phi1 > self.phitoe] = sqrt(2 / self.k * (phi1[phi1 > self.phitoe] - self.C)) - self.D
            h2[:, :] = sqrt(2 / self.k * (phi2 - self.C)) - self.D
        else:
            phi[phi < 0] = 0
            h[phi <= self.phitoe] = sqrt(2 * phi[phi <= self.phitoe] / (self.k * (self.alpha + 1)))
            h[phi > self.phitoe] = sqrt(2 / self.k * (phi[phi > self.phitoe] - self.C)) - self.D
        zeta = -self.D * ones_like(x)
        zeta[phi <= self.phitoe] = -self.alpha * h[phi <= self.phitoe]
        return h, zeta
    
    def plot(self):
        x = hstack((linspace(0, self.d - self.rw, 100), linspace(self.d + self.rw, 3 * self.d, 100)))
        y = linspace(-self.d, self.d, 101)
        x, y = meshgrid(x, y)
        h, zeta = self.head_interface_well(x, y)
        
        # Plot results
        figure(figsize=(8, 8))
        ax1 = subplot(211, aspect='equal')
        contour(x, y, h, 20)
        if self.pumpfresh:
            iwel = count_nonzero(x[0] < self.xmax)
            contour(x[:,:iwel], y[:,:iwel], h[:,:iwel], levels=[self.htoe], colors='k', linewidths=2)
            if h[50, 100] > -self.D:
                title('contour lines of head; head at well: ' + str(round(h[50, 100], 2)) + 'm \n' +
                      'black line is location of interface')
            else:
                title('pump discharge is too large - well pumped dry - solution invalid', color='r')
        else:
            contour(x, y, h, levels=[self.htoe], colors='k', linewidths=2)
            title('pump discharge is too large - well draws salt water - solution invalid', color='r')
        #
        subplot(212, sharex=ax1)
        plot(x[0], h[50])
        if self.pumpfresh:
            itoe = find(zeta[50] == -self.D)[0]
            plot(x[0,:itoe+1], zeta[50,:itoe+1], 'k', lw=2)
        else:
            plot(x[0], zeta[50], 'k', lw=2)
        plot([self.d, self.d], [-self.D, 0], 'k--', lw=3, color=[.8, .8, .8])
        title('cross-section through well')
        show()