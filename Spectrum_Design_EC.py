
import numpy as np
class Spectrum_EC8:

    def __init__(self, ag, S, Tb, Tc, Td, damping = 0.05):
        self.ag = ag
        self.S = S
        self.Tb = Tb
        self.Tc = Tc
        self.Td = Td
        self.damping =damping
        self.damping_factor = np.sqrt(10.0 / (5.0 + 100.0 * self.damping))
        if self.damping_factor < 0.55:
            self.damping_factor = 0.55
        self.get_Sa = np.vectorize(self.get_Sa1)

    def get_Sa1(self, T):
        if T <= self.Tb:
            factor = 1.0 + T / self.Tb * (self.damping_factor * 2.5 - 1.0)
        elif T <= self.Tc:
            factor = self.damping_factor * 2.5
        elif T <= self.Td:
            factor = self.damping_factor * 2.5 * self.Tc / T
        elif T <= 4.0:
            factor = self.damping_factor * 2.5 * self.Tc * self.Td / T**2
        return self.ag * self.S * factor

    def get_Shape(self, num=10):
        Ts1 = np.array([0.0, self.Tb, self.Tc])
        Ts2 = np.linspace(self.Tc, self.Td, num)
        Ts3 = np.linspace(self.Td, 4.0, num)
        Ts = np.hstack((Ts1, Ts2, Ts3))
        return Ts, self.get_Sa(Ts)

