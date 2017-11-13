import numpy as np
import os

from Material_EC import Concrete, Steel_Structural

class CFST:
    #A CFST class according to EC4
    def __init__(self, Materials, Forces, shape, imperfection=True, second_order=True, long_term=False):
        self.shape = shape
        self.long_term = long_term
        self.alpha = 0.21#Imperfection factor
        self.imperfection = imperfection
        self.second_order = second_order
        #Create new materials
        fck, fy, mat_partials = Materials
        if mat_partials==False:
            partial_factor_c, partial_factor_s = 1.5, 1.0
        else:
            partial_factor_c, partial_factor_s = mat_partials
        self.concrete = Concrete(fck=fck, partial_factor=partial_factor_c,long_term=long_term)
        self.steel = Steel_Structural(fy=fy, partial_factor=partial_factor_s)
        #Get the section forces
        self.Ned, self.Med1, self.Med2, self.Ved = Forces
        if self.Ned == 0.0:
            self.Ned = 1.
    
    #Update the member properties
    def update(self):
        #The section geometry properties
        self.Geo_Prop()
        #The section axial resistance
        self.Axial_Resistance()
        #The section shear resistance and steel strength reduction factor
        self.Shear_Resistance()
        self.steel_ratio = self.Aa * self.steel.fyd / self.NplRd#Steel ratio
        #Moment interaction curve
        self.NM_interaction()
        #Axial resistance of a member
        self.Axial_Resistance_Member()
        #Bending resistance of member
        self.Bending_Member()
        #Build a dictionary of member
        self.Build_Dict()
        #Member check
        self.Member_Check()

    #Check the member capacities
    def Member_Check(self):
        Descriptions = []
        #Check the depth to de width ratio
        Descriptions.append('The depth to width ratio is ' + str(self.D2W) +': ') 
        if self.D2W >=0.2 and self.D2W <= 5.0:
            self.is_D2W = True
        else:
            self.is_D2W = False
        Descriptions.append('OK!\n') if self.is_D2W else Descriptions.append('Fail!\n')
        
        #Check the steel contribution ratio:
        Descriptions.append('The steel contribution ratio is ' + str(self.steel_ratio) +': ' )
        if self.steel_ratio >= 0.2 and self.steel_ratio <= 0.9:
            self.is_steel_ratio = True
        else:
            self.is_steel_ratio = False
        Descriptions.append('OK!\n') if self.is_steel_ratio else Descriptions.append('Fail!\n')
        #Check the section slenderness ratio
        Descriptions.append('The section slendernes level is ' + str(self.D2t_level) +': ' )
        if self.shape=='Rectangular':
            if self.D2t_level <=52.:
                self.is_D2t = True
            else:
                self.is_D2t = False

        elif self.shape == 'Circular':
            if self.D2t_level <=90.:
                self.is_D2t = True
            else:
                self.is_D2t = False
        Descriptions.append('OK!\n') if self.is_D2t else Descriptions.append('Fail!\n')
        #Check axial load
        Descriptions.append('Axial load ' + str(self.Nc_member) + ' vs '+ str(self.Ned) +': ' )
        if self.Nc_member >= self.Ned:
            self.is_compression = True
        else:
            self.is_compression = False
        Descriptions.append('OK!\n') if self.is_compression else Descriptions.append('Fail!\n')
        #Check shear load
        Descriptions.append('Shear force ' + str(self.Vpl_a_Rd) + ' vs '+ str(self.Ved) +': ' )
        if self.Vpl_a_Rd >= self.Ved:
            self.is_shear = True
        else:
            self.is_shear = False
        Descriptions.append('OK!\n') if self.is_shear else Descriptions.append('Fail!\n')
        #Check the moment
        Descriptions.append('Moment ' + str(self.M_member) + ' vs '+ str(self.Med_second_order) +': ' )
        if self.M_member >= self.Med_second_order:
            self.is_moment = True
        else:
            self.is_moment = False
        Descriptions.append('OK!\n') if self.is_moment else Descriptions.append('Fail!\n')
        self.output = ''.join(Descriptions)
        

    #The section axial load resistance
    def Axial_Resistance(self):
        self.NplRd = self.Aa * self.steel.fyd + self.Ac * self.concrete.fcd#Section plastic compression capacity
        self.NplRk = self.Aa * self.steel.fy + self.Ac * self.concrete.fck#Section plastic compression capacity calculated with characteristic strength
        self.Nt = self.Aa * self.steel.fyd#Section plastic tensile capacity
        self.Ncr = np.pi**2 * self.EIeff / self.Leff**2#Section buckling critical force
        self.Ncr_eff = np.pi**2 * self.EIeff_2 / self.L**2#Section buckling critical force
        self.L2D_eff = np.sqrt(self.NplRk / self.Ncr)#Relative slenderness
        if self.shape == 'Circular':
            e = max(abs(self.Med1), abs(self.Med2)) / self.Ned
            if e/self.D < 0.1 and self.L2D_eff <= 0.5:
                factor_ao = min(0.25 * (3. + 2. * self.L2D_eff), 1.0)
                factor_co = max(4.9 - 18.5 * self.L2D_eff + 17. * self.L2D_eff**2, 0.0)
                factor_a = factor_ao + (1. - factor_ao) * 10. * e / self.D
                factor_c = factor_co * (1. - 10. * e / self.D)
            else:
                factor_a = 1.0
                factor_c = 0.0
            self.NplRd_uncomfined = self.NplRd
            self.NplRd = factor_a * self.Aa * self.steel.fyd + (1. + factor_c * self.t / self.D * self.steel.fy / self.concrete.fck) * self.Ac * self.concrete.fcd 
        else:
            self.NplRd_uncomfined = self.NplRd


    #The section shear resistance
    def Shear_Resistance(self):
        VaEd = self.Ved#Assume that the shear force is taken by the steel section
        gamma_M0 = 1.0#Partital factor for all section
        self.Vpl_a_Rd = self.Av * self.steel.fy / np.sqrt(3.) / gamma_M0#The plastic shear resistance
        if VaEd >= 0.5 * self.Vpl_a_Rd:
            self.shear_reduction = 1.0 - (2. * VaEd / self.Vpl_a_Rd - 1.0)**2
            self.is_shear_reduce = True
        else:
            self.shear_reduction = 1.0
            self.is_shear_reduce = False

    #Get the moment interaction
    def NM_interaction(self, num=30):
        from scipy.interpolate import interp1d
        fyd = self.steel.fyd * self.shear_reduction
        fcd = self.concrete.fcd
        neutrals = np.linspace(self.neutral_range[0], self.neutral_range[1], num)

        A_sc, cm_sc = self.A_cm_ring(neutrals)
        A_cc, cm_cc = self.A_cm(neutrals)
        A_st, cm_st = self.A_cm_ring(-neutrals)
        A_ct, cm_ct = self.A_cm(-neutrals)
        temp_N = A_sc * fyd + A_cc * fcd - A_st * fyd
        temp_M = A_sc * fyd * cm_sc + A_cc * fcd * cm_cc + A_st * fyd * cm_st
        NM = np.vstack((temp_N, temp_M))
        self.NM = np.hstack(([[-self.Nt], [0.0]], NM, [[self.NplRd], [0.0]]))
        self.f = interp1d(self.NM[0], self.NM[1])
    
    #Get the member axial load capcity accourding to the length
    def Axial_Resistance_Member(self):
        temp_factor = 0.5 * (1. + self.alpha * (self.L2D_eff- 0.2) + self.L2D_eff**2)
        self.axial_reduction = 1. / (temp_factor + np.sqrt(temp_factor**2 - self.L2D_eff**2))
        self.axial_reduction = min(1.0, self.axial_reduction)
        self.Nc_member = self.axial_reduction * self.NplRd

    #Get the member moment capacity with the Ned and Ved
    def Bending_Member(self):
        self.MplRd_N = self.f(self.Ned)#The section mement resistance under Ned, the shear influence was considered in NM interaction
        #Moment reduction
        self.second_order_factors()
        if self.steel.fy <=355.:
            self.alpha_M = 0.9
        else:
            self.alpha_M = 0.8
        self.M_member = self.MplRd_N * self.alpha_M
        #Check if consider imperfection:
        if self.imperfection == True:
            self.imp_level = self.L / 300.
        else:
            self.imp_level = self.imperfection
        #Second order effects
        if self.second_order:
            self.k0 = 1.0 / (1.0 - self.Ned / self.Ncr_eff)
            beta = max(0.44, 0.66 + 0.44 * self.r)
            self.k1 = beta / (1.0 - self.Ned / self.Ncr_eff)
        else:
            self.k0 = 1.0
            self.K1 = 1.0
        #Increase the Med
        Med_max = max(abs(self.Med1), abs(self.Med2))
        Med_second_order = self.k0 * self.Ned * self.imp_level + self.k1 * Med_max 
        self.Med_second_order = max(Med_max, Med_second_order)

    #Get the factors to consider the second moment
    def second_order_factors(self):
        min_M = min(abs(self.Med1), abs(self.Med2))
        max_M = max(abs(self.Med1), abs(self.Med2))
        if max_M == 0.0:
            r = 0.0
        else:
            r = min_M / max_M
        if self.Med1 * self.Med2 >= 0.0:
            self.r = r
        else:
            self.r = -r

    #Build a dictionary of the member properties
    def Build_Dict(self):
        self.Dict = {'Shape':self.shape,
                'Thickness':self.t,
                'Length':self.L,
                'Effective length':self.Leff,
                'Concrete strength':self.concrete.fck,
                'Concrete design strength':self.concrete.fcd,
                'Steel strength':self.steel.fy,
                'Steel design strength':self.steel.fyd,
                'Section compressive capacity':self.NplRd,
                'Section shear capacity':self.Vpl_a_Rd,
                'Section design shear':self.Ved,
                'Axial load-moment interaction':self.NM,
                'Member compressive resistance':self.Nc_member,
                'Design axial load':self.Ned,
                'Section moment resistance':self.MplRd_N,
                'Member moment resistance':self.M_member,
                'Design end moment 1':self.Med1,
                'Design end moment 2':self.Med2,
                'Equivalent design moment':self.Med_second_order,
                'Effective global slenderness':self.L2D_eff,
                'Section depth to width ratio':self.D2W,
                'Section slenderness':self.D2t,
                'Section slenderness level':self.D2t_level,
                'Section steel contribution ratio':self.steel_ratio,
                'Imperfection level':self.imp_level,
                'Long term effect':self.long_term,
                'Second order effect':self.second_order,
                }


class CFST_R(CFST):
    #A rectangular CFST from EC4
    def __init__(self, Geometry, Materials, Forces, imperfection=True, second_order=True, long_term = False):
        self.B, self.H, self.t, self.L, self.Leff = Geometry
        CFST.__init__(self, Materials, Forces, shape = 'Rectangular', imperfection=imperfection, second_order=second_order, long_term=long_term)
        self.update()

    #The section properties
    def Geo_Prop(self):
        self.neutral_range = (self.H/2. - self.t, -self.H/2. + self.t)#neutral axis range for NM interaction curve
        self.A = self.B * self.H#Total section area
        self.Ac = (self.B - 2. * self.t) * (self.H - 2. * self.t)#Concrete area
        self.Aa = self.A - self.Ac#The steel tube area
        self.Av = self.Aa * self.H / (self.B + self.H)#Shear area of rectangular steel tube
        self.Ic = self.Rect_I(self.B - 2. * self.t, self.H - 2. * self.t)#The moment inertia of the concrete
        self.I = self.Rect_I(self.B, self.H)#The moment inertia of the total section
        self.Ia = self.I - self.Ic#The moment inertia of the steel tube
        self.EIeff = self.steel.E * self.Ia + 0.6 * self.concrete.Ecm * self.Ic#The effective section elastic flexural stiffness
        self.EIeff_2 = 0.9 * (self.steel.E * self.Ia + 0.5 * self.concrete.Ecm * self.Ic)#The effective section elastic flexural stiffness consider second order effects
        self.D2W = self.B / self.H#The width to depth ratio
        self.D2t = max(self.B, self.H) / self.t#The cross section slenderness ratio
        self.D2t_level = self.D2t * np.sqrt(self.steel.fy / 235.)#The cross section slenderness level 

    def Rect_I(self,B, H):
        return B * H**3 / 12.

    def A_cm(self, h, B = True, H = True):
        if B==True and H==True:
            B = self.B - 2. * self.t
            H = self.H - 2. * self.t
        A = B * (H / 2. - h)
        cm = H / 2. - (H / 2. - h) / 2.
        return A, cm

    def A_cm_ring(self, h, B = True, H = True, t = True):
        if B==True and H==True and t==True:
            B, H, t = self.B, self.H, self.t
        A1, cm1 = self.A_cm(h, B = B - 2. * t, H = H - 2. * t)
        A, cm = self.A_cm(h, B = B, H = H)
        A2 = A - A1
        cm2 = (A * cm - A1* cm1) / A2
        return A2, cm2


class CFST_C(CFST):
    #A circular CFST from EC4
    def __init__(self, Geometry, Materials, Forces, imperfection=True, second_order=True, long_term = False):
        self.D, self.t, self.L, self.Leff = Geometry
        CFST.__init__(self, Materials, Forces, shape = 'Circular', imperfection=imperfection, second_order=second_order, long_term=long_term)
        self.update()

    #The section properties
    def Geo_Prop(self):
        self.neutral_range = (self.D/2. - self.t*1.01, -self.D/2. + self.t*1.01)#neutral axis range for NM interaction curve
        self.A = np.pi * self.D**2 /4.0#Total section area
        self.Ac = np.pi * (self.D - 2.0 * self.t)**2 /4.0#Concrete area
        self.Aa = self.A - self.Ac#The steel tube area
        self.Av = 2.0 * self.Aa / np.pi#Shear area of rectangular steel tube
        self.Ic = self.Cir_I(self.D - 2.0 * self.t)#The moment inertia of the concrete
        self.I = self.Cir_I(self.D)#The moment inertia of the total section
        self.Ia = self.I - self.Ic#The moment inertia of the steel tube
        self.EIeff = self.steel.E * self.Ia + 0.6 * self.concrete.Ecm * self.Ic#The effective section elastic flexural stiffness
        self.EIeff_2 = 0.9 * (self.steel.E * self.Ia + 0.5 * self.concrete.Ecm * self.Ic)#The effective section elastic flexural stiffness consider second order effects
        self.D2W = 1.0#The width to depth ratio
        self.D2t = self.D / self.t#The cross section slenderness ratio
        self.D2t_level = self.D2t * np.sqrt(self.steel.fy / 235.)#The cross section slenderness level 

    def Cir_I(self, D):
        return np.pi * D**4 / 64.0

    def A_cm(self, h, D = True):
        if D==True:
            D = self.D - 2. * self.t
        theta = np.arccos(2. * h / D)
        A1 = D**2 * theta / 4.
        A2 = D * h * np.sin(theta) / 2.
        A = A1 - A2
        M = (D**2 - 4. * h**2)**1.5 / 12.
        cm = M / A
        return A, cm

    def A_cm_ring(self, h, D = True, t = True):
        if D==True and t==True:
            D, t = self.D, self.t
        A1, cm1 = self.A_cm(h, D = D - 2. * t)
        A, cm = self.A_cm(h, D = D)
        A2 = A - A1
        cm2 = (A * cm - A1* cm1) / A2
        return A2, cm2
