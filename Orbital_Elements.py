from astropy import units as u
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
import numpy as np

class orbital_elements:
    def __init__(self ,
                 inclination ,
                 raan ,
                 eccentricity ,
                 arg_perigee ,
                 mean_anomaly ,
                 a):

        self.inclination = inclination   # [deg]
        self.raan = raan                 # [deg]
        self.eccentricity = eccentricity
        self.arg_perigee = arg_perigee   # [deg]
        self.mean_anomaly = mean_anomaly # [deg]
        self.a = a                       # [m]
        self.mean_anomaly_360 = self.mean_anomaly + 180

    def __repr__(self):
        return (f"{self.__class__.__name__}("
                f"inclination={self.inclination}, "
                f"raan={self.raan}, "
                f"eccentricity={self.eccentricity}, "
                f"arg_perigee={self.arg_perigee}, "
                f"mean_anomaly={self.mean_anomaly}, "
                f"a={self.a})")
    

def orbit_to_elements(x:Orbit):
    # ecc = x.ecc
    # nu_rad = x.nu.to(u.rad).value
    # E = 2 * np.arctan2(np.sqrt((1 - ecc) / (1 + ecc)) * np.tan(nu_rad / 2), 1)
    # E = E.value
    # mean_anomaly = E - ecc * np.sin(E)
    # M_deg = float((mean_anomaly * u.rad).to(u.deg)/u.deg)

    return orbital_elements(inclination=float(x.inc.to(u.deg)/u.deg) ,
                            raan=float(x.raan.to(u.deg)/u.deg) ,
                            eccentricity=float(x.ecc) ,
                            arg_perigee=float(x.argp.to(u.deg)/u.deg) ,
                            # mean_anomaly=M_deg,
                            mean_anomaly=float(x.nu.to(u.deg)/u.deg) ,
                            a=float(x.a.to(u.m)/u.m))

def elements_to_orbit(x:orbital_elements):
    return Orbit.from_classical(attractor=Earth ,
                                a=x.a * u.m ,
                                ecc=x.eccentricity * u.one ,
                                inc=x.inclination * u.deg ,
                                raan=x.raan * u.deg ,
                                argp=x.arg_perigee * u.deg ,
                                nu=(x.mean_anomaly) * u.deg)