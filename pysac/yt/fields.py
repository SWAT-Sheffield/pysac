import numpy as np

from yt.utilities.physical_constants import mu_0

from yt.frontends.gdf.fields import GDFFieldInfo

class SACGDFFieldInfo(GDFFieldInfo):
    known_other_fields = (
        ("velocity_x", ("m/s", ["velocity_x"], None)),
        ("velocity_y", ("m/s", ["velocity_y"], None)),
        ("velocity_z", ("m/s", ["velocity_z"], None)),)
    known_particle_fields = ()

    def __init__(self, *args, **kwargs):
        self._show_field_errors = (("gas" ,"velocity_magnitude"),)
        return super(SACGDFFieldInfo, self).__init__(*args, **kwargs)

    def setup_fluid_fields(self):

        def density(field, data):
            return data['density_pert'] + data['density_bg']
        self.add_field(('gas', 'density'), function=density, units='kg/m**3')

        def mag_field_x(field, data):
            return data['mag_field_x_pert'] + data['mag_field_x_bg']
        self.add_field(('gas','mag_field_x'), function=mag_field_x, units='T')

        def mag_field_y(field, data):
            return data['mag_field_y_pert'] + data['mag_field_y_bg']
        self.add_field(('gas','mag_field_y'), function=mag_field_y, units='T')

        def mag_field_z(field, data):
            return data['mag_field_z_pert'] + data['mag_field_z_bg']
        self.add_field(('gas','mag_field_z'), function=mag_field_z, units='T')

        def internal_energy(field, data):
            return data['internal_energy_pert'] + data['internal_energy_bg']
        self.add_field(('gas','internal_energy'), function=internal_energy, units='Pa')

        def mag_pressure(field, data):
            if data.ds.dimensionality == 2:
                return (data['mag_field_x']**2 +
                        data['mag_field_y']**2) / (2. * mu_0)
            if data.ds.dimensionality == 3:
                return (data['mag_field_x']**2 + data['mag_field_y']**2 +
                        data['mag_field_z']**2) / (2. * mu_0)
        self.add_field(('gas','mag_pressure'), function=mag_pressure, units='Pa',
                       force_override=True)

        def magnetic_field_strength(field, data):
            if data.ds.dimensionality == 2:
                return np.sqrt(data['mag_field_x']**2 +
                        data['mag_field_y']**2)
            if data.ds.dimensionality == 3:
                return np.sqrt(data['mag_field_x']**2 + data['mag_field_y']**2 +
                        data['mag_field_z']**2)
        self.add_field(('gas','magnetic_field_strength'), function=magnetic_field_strength, units='T',
                       force_override=True)

        def thermal_pressure(field, data):
            #p = (\gamma -1) ( e - \rho v^2/2 - B^2/2)
            g1 = data.ds.parameters.get('gamma', 5./3.) -1
            if data.ds.dimensionality == 2:
                kp = (data['density'] * (data['velocity_x']**2 +
                                         data['velocity_y']**2))/2.
            if data.ds.dimensionality == 3:
                kp = (data['density'] * (data['velocity_x']**2 +
                      data['velocity_y']**2 + data['velocity_z']**2))/2.
            return g1 * (data['internal_energy'] - kp - data['mag_pressure'])
        self.add_field(('gas','thermal_pressure'), function=thermal_pressure,
                       units='Pa', force_override=True)

        def alfven_speed(field, data):
            if data.ds.dimensionality == 2:
                r = np.sqrt(data['mag_field_x']**2 +
                            data['mag_field_y']**2) / np.sqrt(mu_0 * data['density'])
            if data.ds.dimensionality == 3:
                r =  np.sqrt(data['mag_field_x']**2 + data['mag_field_y']**2 +
                             data['mag_field_z']**2) / np.sqrt(mu_0 * data['density'])
            return r
        self.add_field(('gas', 'alfven_speed'), function=alfven_speed,
                       units=r'm/s', force_override=True)

        def sound_speed(field, data):
            gamma = data.ds.parameters.get('gamma', 5./3.)
            return np.sqrt((gamma * data['thermal_pressure']) / data['density'])
        self.add_field(('gas','sound_speed'), function=sound_speed,
                       units=r'm/s', force_override=True)

        def plasma_beta(field, data):
            return data['thermal_pressure'] / data['mag_pressure']
        self.add_field(('gas','plasma_beta'), function=plasma_beta,
                       units=r'dimensionless', force_override=True)
