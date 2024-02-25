import numpy as np
from Constants import *
from Extra import clean_3vec

from astropy import units as u
from poliastro.bodies import Earth
from poliastro.twobody import Orbit




class orbital_elements:
    def __init__(self ,
                 inclination ,
                 raan ,
                 eccentricity ,
                 arg_perigee ,
                 mean_anomaly ,
                 a):

        self.inclination = inclination
        self.raan = raan
        self.eccentricity = eccentricity
        self.arg_perigee = arg_perigee
        self.mean_anomaly = mean_anomaly
        self.a = a

    def __repr__(self):
        return (f"{self.__class__.__name__}("
                f"inclination={self.inclination}, "
                f"raan={self.raan}, "
                f"eccentricity={self.eccentricity}, "
                f"arg_perigee={self.arg_perigee}, "
                f"mean_anomaly={self.mean_anomaly}, "
                f"a={self.a})")
    

def orbit_to_elements(x:Orbit):
    return orbital_elements(inclination=float(x.inc.to(u.deg)/u.deg) ,
                            raan=float(x.raan.to(u.deg)/u.deg) ,
                            eccentricity=float(x.ecc) ,
                            arg_perigee=float(x.argp.to(u.deg)/u.deg) ,
                            mean_anomaly=float(x.nu.to(u.deg)/u.deg) ,
                            a=float(x.a.to(u.m)/u.m))

def elements_to_orbit(x:orbital_elements):
    return Orbit.from_classical(attractor=Earth ,
                                a=x.a * u.m ,
                                ecc=x.eccentricity * u.one ,
                                inc=x.inclination * u.deg ,
                                raan=x.raan * u.deg ,
                                argp=x.arg_perigee * u.deg ,
                                nu=x.mean_anomaly * u.deg)


class Satellite:
    def __init__(self , elements):
        self.elements = elements
        self.position , self.velocity = self.find_pos_vel()
        
        self.position = clean_3vec(self.position)
        self.velocity = clean_3vec(self.velocity)

        print(f"\nvel = {self.velocity}")
        print(f"pos = {self.position}\n")
        

    def update(self , dt):
        gravity_acceleration = self.find_grav_acceleration()
        self.velocity += gravity_acceleration * dt
        self.position += self.velocity * dt
        self.elements = self.find_elements()
        

    def find_pos_vel(self):
        orbit = elements_to_orbit(self.elements)
        pos , vel = orbit.rv()
        pos = pos / u.m
        vel = vel / u.m * u.s
        return clean_3vec(pos) , clean_3vec(vel)


    def find_elements(self):
        element = orbit_to_elements(Orbit.from_vectors(Earth , self.position * u.m , self.velocity * u.m / u.s))
        return element


    def find_grav_acceleration(self):
        r = np.linalg.norm(np.array(self.position))
        r_hat = np.array(self.position) / r
        z = self.position[2] # z-component of the position vector

        
        # Gravitational acceleration without J2
        a_magnitude = G * M / (r**2)
        a_gravity = -a_magnitude * r_hat

        # J2 perturbations
        a_r_J2 = -1.5 * J2 * (G * M * R_earth**2) / (r**4) * (1 - 5 * (z/r)**2)
        a_n_J2 = -1.5 * J2 * (G * M * R_earth**2) / (r**4) * (3 * (z/r)**2)

        # Calculate the angular momentum vector
        h_vec = np.cross(self.position, self.velocity)
        h_hat = h_vec / np.linalg.norm(h_vec)

        # Normal component of J2 acceleration should be in the direction of h_hat
        a_n_J2 = -1.5 * J2 * (G * M * R_earth**2) / (r**4) * (3 * (z/r)**2) * h_hat

        # combine the radial and normal components
        a_J2 = a_r_J2 * r_hat + a_n_J2

        # Total gravitational acceleration
        a_total = clean_3vec((a_gravity + a_J2))

        return np.array(a_total)
    
    
    def find_angular_velocity(self):
        """
        Find the angular velocity in [rad s^-1]
        """
        return float(np.sqrt((G*M) / (self.elements.a**3)))











    # def calculate_true_anomaly(self):
    #     # Extract mean anomaly (M) and eccentricity (e) from the orbital elements
    #     M = np.radians(self.elements.mean_anomaly) # M in radians
    #     e = self.elements.eccentricity
        
    #     # Solve Kepler's Equation for Eccentric Anomaly (E) using Newton-Raphson method
    #     E = M # Initial guess for E
    #     tolerance = 1e-10 # Set a tolerance for the solution
    #     while True:
    #         delta_E = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    #         E -= delta_E
    #         if abs(delta_E) < tolerance:
    #             break
        
    #     # Calculate True Anomaly (ν) from Eccentric Anomaly (E)
    #     ν = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
        
    #     # Ensure ν is in the range [0, 2π]
    #     ν = np.mod(ν, 2 * np.pi)
        
    #     return ν


    # def find_velocity_from_elements(self):
        

    #     # Constants
    #     mu = G * M
        
    #     # True anomaly (ν)
    #     ν = self.calculate_true_anomaly()
        
    #     # Perifocal velocity components
    #     v_r = np.sqrt(mu / self.elements.a) * self.elements.eccentricity * np.sin(ν) / (1 + self.elements.eccentricity * np.cos(ν))
    #     v_t = np.sqrt(mu / self.elements.a) * (1 + self.elements.eccentricity * np.cos(ν)) / np.sqrt(1 - self.elements.eccentricity**2)
        
    #     # Rotate to geocentric equatorial frame
    #     Ω = np.radians(self.elements.raan)
    #     i = np.radians(self.elements.inclination)
    #     ω = np.radians(self.elements.arg_perigee)
        
    #     # Rotation matrix from perifocal to geocentric equatorial frame
    #     R = np.array([
    #         [np.cos(Ω) * np.cos(ω) - np.sin(Ω) * np.sin(ω) * np.cos(i), -np.cos(Ω) * np.sin(ω) - np.sin(Ω) * np.cos(ω) * np.cos(i), np.sin(Ω) * np.sin(i)],
    #         [np.sin(Ω) * np.cos(ω) + np.cos(Ω) * np.sin(ω) * np.cos(i), -np.sin(Ω) * np.sin(ω) + np.cos(Ω) * np.cos(ω) * np.cos(i), -np.cos(Ω) * np.sin(i)],
    #         [np.sin(ω) * np.sin(i), np.cos(ω) * np.sin(i), np.cos(i)]
    #     ])
        
    #     # Perifocal velocity vector
    #     v_perifocal = np.array([v_r * np.cos(ν) - v_t * np.sin(ν), v_r * np.sin(ν) + v_t * np.cos(ν), 0])
        
    #     # Velocity in geocentric equatorial frame
    #     velocity = np.dot(R, v_perifocal)
    #     vx, vy, vz = velocity
    #     return [vx, vy, vz]

    # def find_position_from_elements(self):
    #     # Unpack orbital elements
    #     i = np.radians(self.elements.inclination)
    #     Ω = np.radians(self.elements.raan)
    #     e = self.elements.eccentricity
    #     ω = np.radians(self.elements.arg_perigee)
    #     mean_anomaly = np.radians(self.elements.mean_anomaly)
    #     a = self.elements.a
        
    #     # Solve Kepler's Equation for Eccentric Anomaly (E)
    #     E = solve_keplers_equation(mean_anomaly, e)
        
    #     # Calculate True Anomaly (ν)
    #     ν = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
        
    #     # Distance (r)
    #     r = a * (1 - e * np.cos(E))
        
    #     # Position in the orbital plane (perifocal coordinates)
    #     x_p = r * np.cos(ν)
    #     y_p = r * np.sin(ν)
        
    #     # Transform to the geocentric equatorial frame
    #     # Rotation matrix for argument of perigee (ω), inclination (i), and RAAN (Ω)
    #     R = np.array([
    #         [np.cos(Ω) * np.cos(ω) - np.sin(Ω) * np.sin(ω) * np.cos(i), -np.cos(Ω) * np.sin(ω) - np.sin(Ω) * np.cos(ω) * np.cos(i), np.sin(Ω) * np.sin(i)],
    #         [np.sin(Ω) * np.cos(ω) + np.cos(Ω) * np.sin(ω) * np.cos(i), -np.sin(Ω) * np.sin(ω) + np.cos(Ω) * np.cos(ω) * np.cos(i), -np.cos(Ω) * np.sin(i)],
    #         [np.sin(ω) * np.sin(i), np.cos(ω) * np.sin(i), np.cos(i)]
    #     ])
        
    #     position = np.dot(R, np.array([x_p, y_p, 0])) # Position in geocentric equatorial coordinates
    #     return position.tolist()
    

    # def find_elements_from_velocity(self , get_anomaly=False):
    #     """
    #     Calculate orbital elements from state vectors.
    #     Inputs:
    #     - velocity: list or array of velocity components [vx, vy, vz]
    #     - position: list or array of position components [x, y, z]
    #     Outputs:
    #     - Inclination (i)
    #     - Right Ascension of Ascending Node (raan)
    #     - Eccentricity (e)
    #     - Argument of Perigee (omega)
    #     - Mean Anomaly (M0)
    #     - Semi-major axis (a)
    #     """
        
        
    #     r = np.array(self.position)
    #     v = np.array(self.velocity)
    #     mu = G * M
        
    #     h = np.cross(r, v)
    #     N = np.cross([0, 0, 1], h)
    #     e_vec = (np.cross(v, h) / mu) - (r / np.linalg.norm(r))
    #     e = np.linalg.norm(e_vec)
    #     xi = (np.linalg.norm(v)**2 / 2) - (mu / np.linalg.norm(r))
    #     a = -mu / (2 * xi)
    #     i = np.arccos(h[2] / np.linalg.norm(h))
        
    #     # RAAN
    #     if np.linalg.norm(N) > 0:
    #         raan = np.arccos(N[0] / np.linalg.norm(N))
    #         if N[1] < 0:
    #             raan = 2 * np.pi - raan
    #     else:
    #         raan = 0 # 0 by convention for equatorial orbits
        
    #     # Argument of Perigee
    #     if e > 0 and np.linalg.norm(N) > 0:
    #         omega = np.arccos(np.dot(N, e_vec) / (np.linalg.norm(N) * e))
    #         if e_vec[2] < 0:
    #             omega = 2 * np.pi - omega
    #     else:
    #         omega = 0 # 0 for circular orbits
        
    #     # True Anomaly
    #     if e > 0:
    #         v_true_arg = np.clip(np.dot(e_vec, r) / (e * np.linalg.norm(r)), -1, 1)  # Ensure argument is in [-1, 1]
    #         v_true = np.arccos(v_true_arg)*2 - np.pi
    #         if np.dot(r, v) < 0:
    #             v_true = 2 * np.pi - v_true
    #         # Compute mean anomaly for elliptical orbits
    #         E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(v_true / 2))
    #         M0 = E - e * np.sin(E)
    #     else:
    #         v_true = 0
    #         M0 = 0 # For circular orbits, set to 0 or use true anomaly directly
        
    #     # Convert
    #     i = np.degrees(i)
    #     raan = np.degrees(raan)
    #     omega = np.degrees(omega)
    #     M0 = np.degrees(M0) + 180
    #     v_true = np.degrees(v_true) + 180
    #     E = np.degrees(E) + 180
        
    #     if get_anomaly:
    #         return v_true , E , M0
    #     return orbital_elements(inclination=i , raan=raan , eccentricity=e , arg_perigee=omega , mean_anomaly=M0 , a=a)


    
    


        

# def solve_keplers_equation(mean_anomaly, e, tolerance=1e-6):
#         E = mean_anomaly  # Initial guess
#         delta = 1  # Initialize delta to enter the loop
#         while abs(delta) > tolerance:
#             delta = (E - e * np.sin(E) - mean_anomaly) / (1 - e * np.cos(E))
#             E -= delta
#         return E



# # Testing
# if __name__ == "__main__":
#     test_elements = orbital_elements(inclination=0,
#                                              raan=0,
#                                              eccentricity=0,
#                                              arg_perigee=0,
#                                              mean_anomaly=10,
#                                              a=2)
    
#     test_satellite = Satellite(test_elements)

#     get_elements = test_satellite.find_elements_from_velocity()
    
#     print(f"\nInitial elements: {test_elements}")
#     print(f"Final elements  : {get_elements}\n")




