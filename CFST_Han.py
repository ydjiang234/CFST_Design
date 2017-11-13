# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 14:30:50 2016

@author: yjiang
"""

import numpy as np
from CFST_EC4 import CFST_C, CFST_R

class CFST_Han:

    def __init__(self, Material, N, k1=False, reverse_num = 1e3):
        self.fc, self.fy, self.Ec, self.Es = Material
        self.fck = 0.67 * self.fc
        self.N = N
        self.k1 = k1
        self.reverse_num = reverse_num

    def update(self):
        self.section_property()
        self.axial_capacity()
        self.bending()
        self.NM()
        self.pre_backbone()
        self.backbone()

    #NM interaction
    def NM(self, num=30):
        Ns = np.linspace(0., self.Nu)
        Ms = np.array([])
        for N in Ns:
            temp_M = self.find_M(N)
            Ms = np.append(Ms, temp_M)
        self.NM = np.vstack((Ns, Ms))

    #Input the cylic steps, return the response
    def load_with(self, disps):
        cur_unload_path = unload_path((0.,0.), (1.,1.), (2., 2.), (3., 3.))
        output = np.array([0., 0.])
        for i in range(len(disps)):
            cur_d = disps[i]
            if i == 0:
                self.BBP.update_point(cur_d)
            else:
                pre_d = disps[i-1]
                if cur_d == pre_d:
                    1==1
                else:
                    ans = self.point_location(cur_unload_path)
                    if ans == 1:
                        if self.BBP.is_unload(cur_d):
                            P1, P2, P3, P4 = self.BBP.get_unload_path(pre_d)
                            self.BBP.release_point()
                            cur_unload_path = unload_path(P1, P2, P3, P4)
                            temp = cur_unload_path.update_point(cur_d)
                            if temp == False:
                                cur_unload_path.point = False
                                self.BBP.update_point(cur_d)
                        else:
                            self.BBP.update_point(cur_d)
                    elif ans == 2:
                        if cur_unload_path.is_unload(cur_d):
                            temp_unload_d = cur_unload_path.ds[-1]
                            cur_unload_path.release_point()
                            P1, P2, P3, P4 = self.BBP.get_unload_path(temp_unload_d)
                            self.BBP.release_point()
                            cur_unload_path = unload_path(P1, P2, P3, P4)
                            temp = cur_unload_path.update_point(cur_d)
                            if temp == False:
                                cur_unload_path.point = False
                                self.BBP.update_point(cur_d)
                        else:
                            temp = cur_unload_path.update_point(cur_d)
                            if temp == False:
                                cur_unload_path.point = False
                                self.BBP.update_point(cur_d)
            ans = self.point_location(cur_unload_path)
            if ans == 1:
                output = np.vstack((output, [self.BBP.point[0], self.BBP.point[1]]))
            elif ans ==2:
                output = np.vstack((output, [cur_unload_path.point[0], cur_unload_path.point[1]]))
        return output

    def point_location(self, cur_unload_path):
        if cur_unload_path.point == False and self.BBP.point != False:
            return 1
        elif cur_unload_path.point != False and self.BBP.point == False:
            return 2
        else:
            return False

    #Enlarge the protocol
    def enlarge_protocol(self, protocol, num = 1e3):
        from scipy.interpolate import interp1d
        step_no = np.arange(0, len(protocol))
        f = interp1d(step_no, protocol)
        protocols = f(np.linspace(step_no[0], step_no[-1], num))
        return protocols

class CFST_C_Han(CFST_Han):

    def __init__(self, Geometry, Material, N, k1=False, reverse_num=1e3):
        CFST_Han.__init__(self, Material, N, k1=k1, reverse_num=reverse_num)
        self.D, self.t, self.L = Geometry
        self.update()
   
    #Circle area
    def Cir_A(self, D):
        return np.pi * D**2 / 4.

    #Circular moment inertia
    def Cir_I(self, D):
        return np.pi * D**4 / 64.

    #Section property
    def section_property(self):
        D, t, L = self.D, self.t, self.L
        self.A = self.Cir_A(D)
        self.Ac = self.Cir_A(D - 2. * t)
        self.As = self.A - self.Ac
        self.I = self.Cir_I(D)
        self.Ic = self.Cir_I(D - 2. * t)
        self.Is = self.I - self.Ic
        self.lbd = 4. * L / D
        self.alpha = self.As / self.Ac
        self.factor = self.alpha * self.fy / self.fck

    #Axial capacity
    def axial_capacity(self):
        self.fscy = (1.14 + 1.02 * self.factor) * self.fck
        lbdp = np.pi * np.sqrt(self.Es / 0.67 / self.fy)
        lbdo = np.pi * np.sqrt((420. * self.factor + 550.) / self.fscy)
        d = (13000. + 4657. * np.log(235. / self.fy)) * (25. / (self.fck + 5.))**0.3 + (self.alpha / 0.1)**0.05
        e = -d / (lbdp + 35.)**3
        a = (1. + (35. + 2. * lbdp - lbdo) * e) / (lbdp - lbdo)**2
        b = e - 2. * a * lbdp
        c = 1. - a * lbdo**2 - b * lbdo
        if self.lbd <= lbdo:
            self.phi = 1.
        elif self.lbd <= lbdp:
            self.phi = a * self.lbd**2 + b * self.lbd + c
        else:
            self.phi = d / (self.lbd + 35.)**2
        self.Nu = self.phi * self.fscy * self.A
        self.n = self.N / self.Nu

    #Pure bending
    def bending(self):
        rm = 1.1 + 0.48 * np.log(self.factor + 0.1)
        Wscm = np.pi * self.D**3 / 32.
        self.Mu = rm * Wscm * self.fscy
        self.M = self.find_M(self.N)

    #find M when giving a N
    def find_M(self, N):
        so = 1. + 0.18 * self.factor**-1.15
        if self.factor <= 0.4:
            no = 0.5 - 0.2445 * self.factor
        else:
            no = 0.1 + 0.14 * self.factor**-0.84
        Esc = (0.192 * self.fy / 235. + 0.488) * self.fscy / (3.25e-6) / self.fy
        Ne = np.pi**2 * Esc * self.I / self.L**2
        a = 1. - 2. * self.phi**2 * no
        b = (1. - so) / self.phi**3 / no**2
        c = 2. * (so - 1.) / no
        d = 1. - 0.4 * (N / Ne)
        temp = N / self.Nu
        if temp >= 2. * self.phi**3 * no:
            temp_M = self.Mu * d / a * (1. - temp / self.phi)
        else:
            temp_M = self.Mu * d * (1. + b * temp**2 + c * temp)
        return temp_M

    #Prepare for the backbone
    def pre_backbone(self):
        #Ke
        self.Ke = self.Es * self.Is + 0.6 * self.Ec * self.Ic
        #My
        b = self.alpha / 0.1
        c = self.fc / 60.
        if b <= 1.:
            A1 = -0.137
            B1 = -0.468 * b**2 + 0.8 * b + 0.874
            p = 0.566 - 0.789 * b
        else:
            A1 = 0.118 * b - 0.255
            B1 = 1.306 - 0.1 * b
            p = -0.11 * b - 0.113
        if b <= 0.5:
            q = 1.195 - 0.34 * b
        else:
            q = 1.025
        self.My = (A1 * c + B1) / (A1 + B1) / (p * self.n + q) * self.M


    def backbone(self):
        #Ka
        L1 = self.L / 2.
        self.Ka = 3. * self.Ke / L1**3
        if self.k1 != False:
            self.Ka = self.Ka * self.k1 / (self.Ka + self.k1)
        #Py
        if self.n <= 0.3:
            a = 0.96 - 0.002 * self.factor
        else:
            a = (1.4 - 0.34 * self.factor) * self.n + 0.1 * self.factor + 0.54
        if self.factor <= 1.:
            self.Py = a * (0.2 * self.factor + 0.85) / L1 * self.My
        else:
            self.Py = 1.05 * a / L1 * self.My
        self.Pa = 0.6 * self.Py
        #dp
        r = self.lbd / 40.
        s = self.fy / 345.
        if self.n <= 0.5:
            f1 = 1.336 * self.n**2 - 0.044 * self.n + 0.804
        else:
            f1 = 1.126 - 0.02 * self.n
        self.dp = 6.74 * ((np.log(r))**2 - 1.08 * np.log(r) + 3.33) * f1 / (8.7 - s) * self.Py / self.Ka
        #Kd
        if self.n <= 0.7:
            f2 = 3.043 * self.n - 0.21
        else:
            f2 = 0.5 * self.n + 1.57
        if r <= 1.:
            fr = (8. * self.alpha - 8.6) * r + 6. * self.alpha + 0.9
        else:
            fr = (15. * self.alpha - 13.8) * r + 6.1 - self.alpha
        c = self.fc / 60.
        self.Kd = 0.03 * f2 * fr * self.Ka / (c**2 - 3.39 * c + 5.41)
        #Create backbone
        self.BBP = backbone_P(self.Pa, self.Py, self.dp, self.Ka, self.Kd, self.reverse_num, 0.075 * self.L)

    def EC4(self):
        Geometry = self.D, self.t, self.L, self.L
        Force = self.N, 1.0, 0.0, 0.0
        Materials = self.fc, self.fy, False
        temp = CFST_C(Geometry, Materials, Force)
        self.Ke = temp.EIeff
        self.My = temp.f(self.N)
        self.backbone()

class CFST_R_Han(CFST_Han):

    def __init__(self, Geometry, Material, N, k1=False, reverse_num=1e3):
        CFST_Han.__init__(self, Material, N, k1=k1, reverse_num=reverse_num)
        self.B, self.D, self.t, self.L = Geometry
        self.update()

    def Rect_I(self,B, H):
        return B * H**3 / 12.

    #Section property
    def section_property(self):
        B, D, t = self.B, self.D, self.t
        self.beta = D / B
        self.A = B * D
        self.Ac = (B - 2.*t) * (D - 2.*t)
        self.As = self.A - self.Ac
        self.Ic = self.Rect_I(self.B - 2. * self.t, self.D - 2. * self.t)
        self.I = self.Rect_I(self.B, self.D)
        self.Is = self.I - self.Ic
        self.factor = self.As * self.fy / self.Ac / self.fck
        self.alpha = self.As / self.Ac
        self.lbd = 2. * np.sqrt(3.) * self.L / self.D

    #Axial capacity
    def axial_capacity(self):
        self.fscy = (1.18 + 0.85 * self.factor) * self.fck
        lbdo = np.pi * np.sqrt((220. * self.factor + 450.) / self.fscy)
        lbdp = np.pi * np.sqrt(self.Es / 0.62 / self.fy)
        d = (13500. + 4810. * np.log(235. / self.fy)) * (25. / (self.fck + 5.))**0.3 * (self.alpha / 0.1)**0.05
        e = -d / (lbdp + 35.)**3
        a = (1. + (35. + 2. * lbdp - lbdo) * e) / (lbdp - lbdo)**2
        b = e - 2. * a * lbdp
        c = 1. - a * lbdo**2 - b * lbdo
        if self.lbd <= lbdo:
            self.phi = 1.
        elif self.lbd <= lbdp:
            self.phi = a * self.lbd**2 + b * self.lbd + c
        else:
            self.phi = d / (self.lbd + 35.)**2
        self.Nu = self.phi * self.fscy * self.A
        self.n = self.N / self.Nu

    #Pure bending
    def bending(self):
        B, D, t = self.B, self.D, self.t
        gamma_m = 1.04 + 0.48 * np.log(self.factor + 0.1)
        Wscm = B * D**2 / 6.
        self.Mu = gamma_m * Wscm * self.fscy
        self.M = self.find_M(self.N)

    #find M when giving a N
    def find_M(self, N):
        Ne = np.pi**2 * (self.Es*self.Is + self.Ec * self.Ic) / self.L**2
        if self.factor <= 0.4:
            eta_o = 0.5 - 0.3175 * self.factor
        else:
            eta_o = 0.1 + 0.13 * self.factor**-0.81
        s_o = 1. + 0.14 * self.factor**-1.3

        a = 1. - 2. * self.phi**2 * eta_o
        b = (1. - s_o)/ (self.phi**3 * eta_o**2)
        c = 2. * (s_o - 1.) / eta_o
        d = 1. - 0.25 * N / Ne
        limit = 2. * self.phi**3 * eta_o
        x = N / self.Nu
        if x >= limit:
            M = self.Mu * d / a * (1. - x / self.phi)
        else:
            M = self.Mu * d * (1.0 + b * x**2 + c * x)
        return M

    #Prepare for the backbone
    def pre_backbone(self):
        #Ke
        self.Ke = (self.Es * self.Is + 0.2 * self.Ec * self.Ic)
        #My
        self.My = self.M


    def backbone(self):
        #Ka
        self.Ka = 24. * self.Ke / self.L**3
        if self.k1 != False:
            self.Ka = self.Ka * self.k1 / (self.Ka + self.k1)
        #Py
        if self.n <= 0.4:
            self.Py = (2.5 * self.n**2 - 0.75 * self.n + 1.) * self.M / self.L * 2.
        else:
            self.Py = (0.63 * self.n + 0.848) * self.M / self.L * 2.
        #dp
        self.dp = (1.7 + self.n + 0.5 * self.factor) *self.Py / self.Ka
        #Kd
        self.Kd = -9.83 * self.n**1.2 * self.lbd**0.75 * self.fy / self.Es / self.factor * self.Ka
        #Pa
        self.Pa = 0.6 * self.Py
        self.da = self.Pa / self.Ka
        self.BBP = backbone_P(self.Pa, self.Py, self.dp, self.Ka, self.Kd, self.reverse_num, self.L * 0.075)

    def EC4(self):
        Geometry = self.B, self.D, self.t, self.L, self.L
        Force = self.N, 1.0, 0.0, 0.0
        Materials = self.fc, self.fy, False
        temp = CFST_R(Geometry, Materials, Force)
        self.Ke = temp.EIeff
        self.My = temp.f(self.N)
        self.backbone()

class backbone_P:

    def __init__(self,Pa,Pb,db,Ka,Kd, reverse_num, limit):
        from scipy.interpolate import interp1d
        P1 = (Pa/Ka, Pa)
        P2 = (db, Pb)
        if db + Pb / abs(Kd) > limit:
            P3 = (limit, (limit - db) * Kd + Pb)
        else:
            P3 = (db + Pb / abs(Kd), 0.)
        self.Ka = Ka
        self.ds = np.array([-P3[0], -P2[0], -P1[0], P1[0], P2[0], P3[0]])
        self.Ps = np.array([-P3[1], -P2[1], -P1[1], P1[1], P2[1], P3[1]])
        self.f = interp1d(self.ds, self.Ps)
        self.Prs = -0.2 * self.Ps
        self.fr = interp1d(self.ds, self.Prs)
        self.drs = np.linspace(self.ds[0], self.ds[-1], reverse_num)
        self.Prs = self.fr(self.drs)
        self.point = False
        self.linear_range = (-P1[0], P1[0])

    #Is the point in linear range
    def is_linear(self):
        if self.point != False:
            if self.point[0]<=self.linear_range[1] and self.point[0]>= self.linear_range[0]:
                return True
            else:
                return False

    #Is unload in next step
    def is_unload(self, cur_d):
        if self.point != False:
            pre_d = self.point[0]
            if not self.is_linear():
                if (cur_d > pre_d and pre_d>0) or (cur_d<pre_d and pre_d < 0.):
                    return False
                else:
                    return True
            else:
                return False

    #Get unload path
    def get_unload_path(self,pre_d):
        pre_P = self.f(pre_d)
        P1 = (pre_d, pre_P)
        P4 = (-pre_d, -pre_P)
        P3 = (-pre_P/self.Ka, -pre_P)
        Ks = abs((pre_P - self.Prs) / (pre_d - self.drs) - self.Ka)
        min_id = np.argmin(Ks)
        P2 = (self.drs[min_id], self.Prs[min_id])
        return P1, P2, P3, P4

    #Put the point to the Backbone
    def update_point(self, cur_d):
        cur_P = self.f(cur_d)
        self.point = (cur_d, cur_P)

    def release_point(self):
        self.point = False

class unload_path:

    def __init__(self, P1, P2, P3, P4):
        from scipy.interpolate import interp1d
        self.ds = np.array([P1[0], P2[0], P3[0], P4[0]])
        self.Ps = np.array([P1[1], P2[1], P3[1], P4[1]])
        if self.ds[0] < self.ds[-1]:
            self.f = interp1d(self.ds, self.Ps)
            self.direction = -1
        else:
            self.f = interp1d(self.ds[::-1], self.Ps[::-1])
            self.direction = 1
        self.point = False
        self.linear_range = (min(P1[0], P3[0]), max(P1[0], P3[0]))
        self.path_range = (min(self.ds), max(self.ds))
    
    def update_point(self, cur_d):
        if cur_d < self.path_range[0] or cur_d > self.path_range[1]:
            return False
        else:
            self.point = (cur_d, self.f(cur_d))
            return True

    def release_point(self):
        self.point = False

    def is_linear(self):
        if self.point != False:
            pre_d = self.point[0]
            if pre_d <= self.linear_range[1] and pre_d >= self.linear_range[0]:
                return True
            else:
                return False
    def is_unload(self, cur_d):
        if self.point != False:
            pre_d = self.point[0]
            if not self.is_linear():
                if cur_d > pre_d and self.direction == 1:
                    return True
                elif cur_d < pre_d and self.direction == -1:
                    return True
                else:
                    return False
            else:
                return False

