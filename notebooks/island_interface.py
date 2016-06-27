from __future__ import division
from __future__ import print_function
from pylab import *

class IslandInterface:
    def __init__(self, k=10, D=50, c=0, rhof=1000, rhos=1025, L=1000, N=0.001, Nstreamlines=None):
        self.k = k
        self.D = D
        self.c = c
        self.rhof = rhof
        self.rhos = rhos
        self.L = L
        self.N = N
        self.Nstreamlines = Nstreamlines
        #
        self.alpha = self.rhof / (self.rhos - self.rhof)
        self.phitoe = 0.5 * self.k * (self.alpha + 1) / (self.alpha ** 2) * self.D ** 2
        self.C = -0.5 * self.k * (self.alpha + 1) / self.alpha * self.D ** 2
        #
        if self.c > 0:
            self.lab = sqrt(self.k * self.D * self.c)
            self.nu = (self.rhos - self.rhof) / self.rhof
            self.grad = self.N * self.L / (self.k * self.D)
            self.mu =  self.grad * self.lab / self.D / self.nu
            self.h0 = self.nu * self.D * (3 * self.mu ** 2 / 2) ** (1 / 3)
            self.Loutflow = (18 * self.mu) ** (1 / 3) * self.lab
    def plot(self):
        if self.c == 0:
            x = linspace(-self.L, self.L, 100)
            phi = -self.N / 2 * (x ** 2 - self.L ** 2)
            h =  nan * ones(len(x))
            h[phi <= self.phitoe] = sqrt(2 * phi[phi <= self.phitoe] / (self.k * (self.alpha + 1)))
            h[phi > self.phitoe] = sqrt(2 / self.k * (phi[phi > self.phitoe] - self.C)) - self.D
            zeta = -self.D * ones(len(x))
            zeta[phi <= self.phitoe] = -self.alpha * h[phi <= self.phitoe]
            # Plot results
            figure(figsize=(10, 5))
            plot(x, h, 'b')
            plot(x, zeta, 'r')
            # Plot streamlines
            xs = np.vstack((x, x))
            ys = np.vstack((h, zeta))
            Qs = np.vstack((self.N * x, zeros(len(x))))
            if self.Nstreamlines is not None and self.Nstreamlines > 0: 
                contour(xs, ys, Qs, self.Nstreamlines)
        else:
            x = linspace(-self.L, self.L, 100)
            phicoast = 0.5 * self.k * (self.alpha + 1) * self.h0 ** 2
            phi = -self.N / 2 * (x ** 2 - self.L ** 2) + phicoast
            h =  nan * ones(len(x))
            h[phi <= self.phitoe] = sqrt(2 * phi[phi <= self.phitoe] / (self.k * (self.alpha + 1)))
            h[phi > self.phitoe] = sqrt(2 / self.k * (phi[phi > self.phitoe] - self.C)) - self.D
            zeta = -self.D * ones(len(x))
            zeta[phi <= self.phitoe] = -self.alpha * h[phi <= self.phitoe]
            # Outflow zone
            x2 = linspace(0, self.Loutflow, 100)
            phi2 = (x2 / self.lab - self.Loutflow / self.lab) ** 2 / 6
            h2 = self.nu * self.D * phi2
            zeta2 = -h2 / self.nu
            # Plot results
            figure(figsize=(10, 5))
            plot(x, h, 'b')
            plot(x2 + self.L, h2, 'b')
            plot(-x2 - self.L, h2, 'b')
            plot(x, zeta, 'r')
            plot(x2 + self.L, zeta2, 'r')
            plot(-x2 - self.L, zeta2, 'r')
            # Plot streamlines
            xh = np.hstack(((-x2 - self.L)[::-1], x, x2 + self.L))
            Qoutflow = -self.k * (0 - zeta2) * self.D * self.nu * (x2 / self.lab - self.Loutflow / self.lab) / (3 * self.lab)
            Qh = np.hstack((-Qoutflow[::-1], self.N * x, Qoutflow))
            htop = np.hstack((h2[::-1], h, h2))
            zetabot = np.hstack((zeta2[::-1], zeta, zeta2))
            xs = np.vstack((xh, xh))
            ys = np.vstack((htop, zetabot))
            Qs = np.vstack((Qh, zeros(len(xh))))
            if self.Nstreamlines is not None and self.Nstreamlines > 0:
                contour(xs, ys, Qs, self.Nstreamlines)
        show()
        
        
#island = IslandInterface(k=10, D=1000, c=0, rhof=1000, rhos=1025, L=1000, N=0.001, Nstreamlines=20)
#island.plot()