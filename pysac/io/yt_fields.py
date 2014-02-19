# -*- coding: utf-8 -*-
import yt.mods as yt
import numpy as np

mu0 = 1.25663706e-6
gamma = 1.6666

@yt.derived_field(take_log=False, units=r'kg m^{-3}')
def density(field, data):
    return data['density_pert'] + data['density_bg']

@yt.derived_field(take_log=False, units=r'T')
def mag_field_x(field, data):
    return data['mag_field_x_pert'] + data['mag_field_x_bg']

@yt.derived_field(take_log=False, units=r'T')
def mag_field_y(field, data):
    return data['mag_field_y_pert'] + data['mag_field_y_bg']

@yt.derived_field(take_log=False, units=r'T')
def mag_field_z(field, data):
    return data['mag_field_z_pert'] + data['mag_field_z_bg']

@yt.derived_field(take_log=False, units=r'T')
def mag_field_magnitude(field, data):
    return np.sqrt(data['mag_field_x']**2 + data['mag_field_y']**2 + data['mag_field_z']**2)

@yt.derived_field(take_log=False, units=r'T')
def mag_field_pert_magnitude(field, data):
    return np.sqrt(data['mag_field_x_pert']**2 + data['mag_field_y_pert']**2 +
                    data['mag_field_z_pert']**2)

@yt.derived_field(take_log=False, units=r'T')
def velocity_magnitude(field, data):
    return np.sqrt(data['velocity_x']**2 + data['velocity_x']**2 +
                    data['velocity_x']**2)

@yt.derived_field(take_log=False, units=r'Pa')
def internal_energy(field, data):
    return data['internal_energy_pert'] + data['internal_energy_bg']

@yt.derived_field(take_log=False, units=r'Pa')
def mag_pressure(field, data):
    return (data['mag_field_x']**2 + data['mag_field_y']**2 + data['mag_field_z']**2) / (2. * mu0)

@yt.derived_field(take_log=False, units=r'Pa')
def thermal_pressure(field, data):
    #p = (\gamma -1) ( e - \rho v^2/2 - B^2/2)
    g1 = gamma -1 #(header['eqpar'][0]-1)
    kp = (data['density'] * (data['velocity_x']**2 + data['velocity_y']**2 + data['velocity_z']**2))/2.
    return g1 * (data['internal_energy'] - kp - data['mag_pressure'])

@yt.derived_field(take_log=False, units=r'm s^{-1}')
def alfven_speed(field, data):
    return np.sqrt(data['mag_field_x']**2 + data['mag_field_y']**2 + data['mag_field_z']**2) / np.sqrt(data['density'])

@yt.derived_field(take_log=False, units=r'm s^{-1}')
def sound_speed(field, data):  
    return np.sqrt((gamma * data['thermal_pressure']) / data['density'])

@yt.derived_field(take_log=False, units=r'')
def plasma_beta(field, data):  
    return data['mag_pressure'] / data['thermal_pressure']
#    
#@yt.derived_field(take_log=False, units=r'Pa')
#def wave_flux_x(field, data):
#    Bb = np.array([f.w_dict['bg3'], f.w_dict['bg2'], f.w_dict['bg1']])
#    Bp = np.array([f.w_dict['b3'], f.w_dict['b2'], f.w_dict['b1']])
#    V = np.array([f.w_sac['v3'], f.w_sac['v2'], f.w_sac['v1']])
#    
#    #Calculate wave flux
#    Fp = 0.25*np.pi * (np.sum(Bb*Bp, axis=0)[None] * V) - (np.sum(V*Bp, axis=0)[None] * Bb)
#    Fa = pk[None]*V
#    
#    Fwave = Fa + Fp
    
#    def get_total_p(self):
#        if self.header['ndim'] == 3:
#           gamma = self.header['eqpar'][0]
#           
#           vtot2 = (self.w_sac['v1']**2 + self.w_sac['v2']**2 + self.w_sac['v3']**2)
#           therm = self.w[self.w_["e"]] - (self.w_sac["rho"] * vtot2) / 2.
#           
#           Bpert = self.w[self.w_['b1']] + self.w[self.w_['b2']] + self.w[self.w_['b3']]
#           Bpert2 = self.w[self.w_['b1']]**2 + self.w[self.w_['b2']]**2 + self.w[self.w_['b3']]**2
#           Bback = self.w[self.w_['bg1']] + self.w[self.w_['bg2']] + self.w[self.w_['bg3']]
#           mag = Bback * Bpert + (Bpert2 / 2.)
#           
#           return (gamma - 1) * therm - (gamma - 2) * mag
#        else:
#            raise NotImplementedError("This Dosen't work for 2D yet, go fix")
#    
#    def get_temp(self,p=None):
#        if not(p):
#            p = self.get_thermalp()
#        T = (p * 1.2) / (8.3e3 * self.w_sac['rho'])
#        return T
#    
#    def get_bgtemp(self):
#        print "WARNING: Background Temprature will not work if inital conditions are not V=0"
#        if self.header['ndim'] == 3:
#            kp = 0.0#(self.w[self.w_["rhob"]] * (self.w_sac['v1']**2 + self.w_sac['v2']**2 + self.w_sac['v3']**2))/2.
#            mp = (self.w[self.w_["bg1"]]**2 + self.w[self.w_["bg2"]]**2 + self.w[self.w_["bg3"]]**2) / 2.
#            T = self.w[self.w_["eb"]] - kp - mp
#        else:
#            kp = 0.0#(self.w[self.w_["rhob"]] * (self.w_sac['v1']**2 + self.w_sac['v2']**2))/2.
#            mp = (self.w[self.w_["bg1"]]**2 + self.w[self.w_["bg2"]]**2) / 2.
#            T = self.w[self.w_["eb"]] - kp - mp
#        return T
#    
#    def get_va(self):
#        return (np.sqrt(self.w_sac['b1']**2 + self.w_sac['b2']**2
#                        + self.w_sac['b3']**2) / np.sqrt(self.w_sac['rho']))
#        #return (abs(self.w_sac['b1']) + abs(self.w_sac['b2']) + abs(self.w_sac['b3'])) / sqrt(self.w_sac['rho'])
#    
#    def get_cs(self,p=None):
#        if not p:
#            p = self.get_thermalp()
#        g1 = self.header['eqpar'][0]
#        return np.sqrt((g1 * p) / self.w_sac['rho'])