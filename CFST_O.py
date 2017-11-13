# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 2017

@author: yjiang
"""

class CFST_O:

    def __init__(self, D, t, fc, fy):
        self.D = D
        self.t = t
        self.fc = fc
        self.fy = fy
        self.update()

    def update(self):
        import numpy as np
        from scipy.interpolate import interp1d
        self.A = self.D**2 * np.sin(np.pi/4.)
        self.Ac = (self.D - 2. * self.t)**2 * np.sin(np.pi/4.)
        self.As = self.A - self.Ac
        self.Nc = self.As * self.fy + self.Ac * self.fc
        self.Nt = self.As * self.fy
        B1 = self.D * np.sin(3. * np.pi / 8.)
        B2 = self.D * np.cos(3. * np.pi / 8.)
        h = B2 / 2.**0.5
        I1 = B1 * B2**2 / 12.
        A2, I2, cm2 = self.trapezium(B1, B2, h)
        I3 = I2 + A2 * B2**2 / 4.
        self.I = I1 + 2. * I3
        self.NM = self.N_M_interactive()
        self.f = interp1d(self.NM[:,0][::-1], self.NM[:,1][::-1])

    def trapezium(self, B1, B2, h):
        A2 = (B1 + B2) * h / 2.
        I2 = B1 * h**3 * (1. + 3. * B2 / B1)
        H = B1 * h / (B1 - B2)
        h1 = H - h
        A = 0.5 * B1 * H
        A1 = 0.5 * B2 * h1
        cm = H / 3.
        cm1 = h + h1 / 3.
        cm2 = (A * cm - A1 * cm1) / A2
        return A2, I2, cm2
        
    def A_cm_Octa(self, D, h):
        import numpy as np
        B1 = D * np.sin(3. * np.pi / 8.)
        B2 = D * np.cos(3. * np.pi / 8.)
        h1 = B2 / 2.**0.5
        if h <= 0.5 * B2 and h >= -0.5 * B2:
            A1, I1, cm1 = self.trapezium(B1, B2, h1)
            cm1 = cm1 + 0.5 * B2
            A2 = (0.5 * B2 - h) * B1
            cm2 = 0.5 * B2 - 0.5 * (0.5 * B2 - h)
            A = A1 + A2
            cm = (A1 * cm1 + A2 * cm2) / A
        elif h > 0.5 * B2:
            temp_B2 = B2
            temp_h = 0.5 * B1 - h
            temp_B1 = temp_B2 + 2. * temp_h
            A, I ,cm = self.trapezium(temp_B1, temp_B2, temp_h)
            cm = cm + h
        elif h < -0.5 * B2:
            A1, I1, cm1 = self.trapezium(B1, B2, h1)
            cm1 = cm1 + 0.5 * B2
            A2 = B1 * B2
            cm2 = 0.
            temp_B1 = B1
            temp_h = -h - 0.5 * B2
            temp_B2 = temp_B1 - 2. * temp_h
            A3, I3 ,cm3 = self.trapezium(temp_B1, temp_B2, temp_h)
            cm3 = cm3 + 0.5 * B2
            A = A1 + A2 + A3
            cm = (A1 * cm1 + A2 * cm2 - A3 * cm3) / A
        return A, cm

    def A_cm_Octa_ring(self, D, t, h):
        A1, cm1 = self.A_cm_Octa(D - 2. * t, h)
        A, cm = self.A_cm_Octa(D, h)
        A2 = A - A1
        cm2 = (A * cm - A1* cm1) / A2
        return A2, cm2
        
    def N_M_interactive(self, num = 30):
        import numpy as np
        B1 = self.D * np.sin(3. * np.pi / 8.)
        B2 = self.D * np.cos(3. * np.pi / 8.)
        neutrals = np.linspace(-B1/2.0 + self.t, B1/2.0 - self.t, num)
        neutrals = np.append(neutrals, [-0.5 * B2, 0.5 * B2])
        neutrals = np.sort(neutrals, axis=None)
        NM = np.array([self.Nc, 0.0])     
        for neutral in neutrals:
            A_sc, cm_sc = self.A_cm_Octa_ring(self.D, self.t, neutral)
            A_cc, cm_cc = self.A_cm_Octa(self.D-2.0*self.t, neutral)
            A_st, cm_st = self.A_cm_Octa_ring(self.D, self.t, -neutral)
            A_ct, cm_ct = self.A_cm_Octa(self.D-2.0*self.t, -neutral)
            temp_N = A_sc * self.fy + A_cc * self.fc - A_st * self.fy
            temp_M = A_sc * self.fy * cm_sc + A_cc * self.fc * cm_cc + A_st * self.fy * cm_st
            NM = np.vstack((NM, [temp_N, temp_M]))
        NM = np.vstack((NM, [-self.Nt, 0.0]))
        return NM
    
    def find_M(self, N):
        return self.f(N)
    
    def find_M_p(self, per):
        return self.f(self.Nc*per)
