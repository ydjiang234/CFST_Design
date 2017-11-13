#Contents material class in Eurocode

import os
import numpy as np

class Concrete:
#A concrete class from EC2
    def __init__(self, fck, fck_cube=0.0, partial_factor=1.5, alpha_cc=1.0, long_term=False,label='Concrete'):
        self.fck = fck#Characteristic compressive cylinder strength of concrete at 28 days
        self.fck_cube = fck_cube #Concrete characteristic cube strength
        self.partial_factor = partial_factor #Partial factor for concrete
        self.alpha_cc = alpha_cc#Long term reduction
        self.label = label#Description for the concrete
        self.update()
    
    def update(self):
        self.fcm = self.fck + 8.#Mean value of concrete cylinder compressive strength
        self.Ecm = 22. * (self.fcm / 10.)**0.3 * 1000.#Secant modulus of elasticity of concrete
        self.fcd = self.fck / self.partial_factor#Design value of concrete compressive strength



class Steel_Structural:
#A structural steel class from EC3
    def __init__(self, fy, partial_factor = 1.0, label='Structual steel'):
        self.fy = fy#Characteristic yield strength
        self.partial_factor = partial_factor#Partial factor for structural steel
        self.label = label
        self.update()
    
    def update(self):
        self.E = 210.e3#modulus of elasticity
        self.poisson = 0.3#Poisson ratio in elastic stage
        self.G = self.E / 2. / (1. + self.poisson)#shear modulus
        self.fyd = self.fy / self.partial_factor#Design steel strength
