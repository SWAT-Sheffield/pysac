import numpy as np

from yt.frontends.gdf.fields import GDFFieldInfo

class SACGDFFieldInfo(GDFFieldInfo):
    known_other_fields = (
        ("density", ("g/cm**3", ["density"], None)),
        ("velocity_x", ("cm/s", ["velocity_x"], None)),
        ("velocity_y", ("cm/s", ["velocity_y"], None)),
        ("velocity_z", ("cm/s", ["velocity_z"], None)),
        ("mag_field_x", ("gauss", ["mag_field_x"], None)),
        ("mag_field_y", ("gauss", ["mag_field_y"], None)),
        ("mag_field_z", ("gauss", ["mag_field_z"], None)),)
    known_particle_fields = ()

    def __init__(self, *args, **kwargs):
        self._show_field_errors = (("gas" ,"velocity_magnitude"),)
        return super(SACGDFFieldInfo, self).__init__(*args, **kwargs)

    def setup_fluid_fields(self):
        mu0 = np.pi * 4

        def density(field, data):
            return data['density_pert'] + data['density_bg']
        self.add_field(('gdf', 'density'), function=density)

        def mag_field_x(field, data):
            return data['mag_field_x_pert'] + data['mag_field_x_bg']
        self.add_field(('gdf','mag_field_x'), function=mag_field_x)

        def mag_field_y(field, data):
            return data['mag_field_y_pert'] + data['mag_field_y_bg']
        self.add_field(('gdf','mag_field_y'), function=mag_field_y)

        def mag_field_z(field, data):
            return data['mag_field_z_pert'] + data['mag_field_z_bg']
        self.add_field(('gdf','mag_field_z'), function=mag_field_z)

        def internal_energy(field, data):
            return data['internal_energy_pert'] + data['internal_energy_bg']
        self.add_field(('gdf','internal_energy'), function=internal_energy)

        def mag_pressure(field, data):
            if data.ds.dimensionality == 2:
                return (data['mag_field_x']**2 +
                        data['mag_field_y']**2) / (2. * mu0)
            if data.ds.dimensionality == 3:
                return (data['mag_field_x']**2 + data['mag_field_y']**2 +
                        data['mag_field_z']**2) / (2. * mu0)
        self.add_field(('gdf','mag_pressure'), function=mag_pressure)

        def thermal_pressure(field, data):
            #p = (\gamma -1) ( e - \rho v^2/2 - B^2/2)
            g1 = data.ds.parameters.get('gamma',[1.66666])[0] -1
            if data.ds.dimensionality == 2:
                kp = (data['density'] * (data['velocity_x']**2 +
                                         data['velocity_y']**2))/2.
            if data.ds.dimensionality == 3:
                kp = (data['density'] * (data['velocity_x']**2 +
                      data['velocity_y']**2 + data['velocity_z']**2))/2.
            return g1 * (data['internal_energy'] - kp - data['mag_pressure'])
        self.add_field(('gdf','thermal_pressure'), function=thermal_pressure)

        def alfven_speed(field, data):
            if data.ds.dimensionality == 2:
                return np.sqrt(data['mag_field_x']**2 +
                               data['mag_field_y']**2) / np.sqrt(mu0 * data['density'])
            if data.ds.dimensionality == 3:
                return np.sqrt(data['mag_field_x']**2 + data['mag_field_y']**2 +
                               data['mag_field_z']**2) / np.sqrt(mu0 * data['density'])
        self.add_field(('gdf', 'alfven_speed'), function=alfven_speed)

        def sound_speed(field, data):
            gamma = data.ds.parameters.get('gamma',[1.66666])[0]
            return np.sqrt((gamma * data['thermal_pressure']) / data['density'])
        self.add_field(('gdf','sound_speed'), function=sound_speed)

        def plasma_beta(field, data):
            return data['mag_pressure'] / data['thermal_pressure']
        self.add_field(('gdf','plasma_beta'), function=plasma_beta, units=r'')
