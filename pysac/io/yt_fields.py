"""
A set of dervied fields for yt 2.x which combine pertubation and background components
and define magnitudes, as well as calculate things like characteristic speeds.

Note: These use yt 3.x like field naming conventions
"""
import warnings

import numpy as np

import yt

if yt.__version__.startswith('2'):

    import yt.mods as yt
    __all__ = ['density', 'velocity_magnitude', 'internal_energy',
            'mag_pressure', 'thermal_pressure', 'alfven_speed', 'sound_speed',
            'plasma_beta', 'mag_field_x', 'mag_field_y', 'mag_field_z',
            'mag_field_magnitude', 'mag_field_pert_magnitude']

    #mu0 = 1.25663706e-6
    mu0 = np.pi * 4
    gamma = 5./3.

    @yt.derived_field(take_log=False, units=r'g cm^{-3}')
    def density(field, data):
        return data['density_pert'] + data['density_bg']

    @yt.derived_field(take_log=False, units=r'G')
    def mag_field_x(field, data):
        return data['mag_field_x_pert'] + data['mag_field_x_bg']

    @yt.derived_field(take_log=False, units=r'G')
    def mag_field_y(field, data):
        return data['mag_field_y_pert'] + data['mag_field_y_bg']

    @yt.derived_field(take_log=False, units=r'G')
    def mag_field_z(field, data):
        return data['mag_field_z_pert'] + data['mag_field_z_bg']

    @yt.derived_field(take_log=False, units=r'G')
    def mag_field_magnitude(field, data):
        return np.sqrt(data['mag_field_x']**2 + data['mag_field_y']**2 + data['mag_field_z']**2)

    @yt.derived_field(take_log=False, units=r'G')
    def mag_field_pert_magnitude(field, data):
        return np.sqrt(data['mag_field_x_pert']**2 + data['mag_field_y_pert']**2 +
                        data['mag_field_z_pert']**2)

    @yt.derived_field(take_log=False, units=r'cm s^{-1}')
    def velocity_magnitude(field, data):
        return np.sqrt(data['velocity_x']**2 + data['velocity_y']**2 +
                    data['velocity_z']**2)

    @yt.derived_field(take_log=False, units=r'Ba')
    def internal_energy(field, data):
        return data['internal_energy_pert'] + data['internal_energy_bg']

    @yt.derived_field(take_log=False, units=r'Ba')
    def mag_pressure(field, data):
        return (data['mag_field_x']**2 + data['mag_field_y']**2 + data['mag_field_z']**2) / (2. * mu0)

    @yt.derived_field(take_log=False, units=r'Ba')
    def thermal_pressure(field, data):
        #p = (\gamma -1) ( e - \rho v^2/2 - B^2/2)
        g1 = gamma -1 #(header['eqpar'][0]-1)
        kp = (data['density'] * (data['velocity_x']**2 + data['velocity_y']**2 + data['velocity_z']**2))/2.
        return g1 * (data['internal_energy'] - kp - data['mag_pressure'])

    @yt.derived_field(take_log=False, units=r'cm s^{-1}')
    def alfven_speed(field, data):
        return np.sqrt(data['mag_field_x']**2 + data['mag_field_y']**2 + data['mag_field_z']**2) / np.sqrt( mu0 * data['density'])

    @yt.derived_field(take_log=False, units=r'cm s^{-1}')
    def sound_speed(field, data):
        return np.sqrt((gamma * data['thermal_pressure']) / data['density'])

    @yt.derived_field(take_log=False, units=r'')
    def plasma_beta(field, data):
        return data['mag_pressure'] / data['thermal_pressure']

    @yt.derived_field(take_log=False, units=r'K')
    def temperature(field, data):
        warnings.warn("Hard coded numbers", Warning)
        return (data['thermal_pressure'] * 1.2) / (8.3e3 * data['density'])
