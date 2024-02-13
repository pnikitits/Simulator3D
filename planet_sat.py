import numpy as np
from Extra import G , M , normalize_vector

J2 = 0#1.08263e-3 # 0.2  # Earth's second gravitational harmonic coefficient

class Satellite:
    def __init__(self , elements):
        self.elements = elements # satelitte's orbital elements
        self.position = self.find_position_from_elements() # position : [x,y,z]
        self.velocity = self.find_velocity_from_elements() # velocity : [vx,vy,vz]


    def update(self , dt):
        gravity_acceleration = self.find_grav_acceleration()
        self.velocity += gravity_acceleration * dt
        self.position += self.velocity * dt
        self.elements = self.find_elements_from_velocity()
        

    def calculate_true_anomaly(self):
        # Extract mean anomaly (M) and eccentricity (e) from the orbital elements
        M = np.radians(self.elements.mean_anomaly)  # Ensure M is in radians
        e = self.elements.eccentricity
        
        # Solve Kepler's Equation for Eccentric Anomaly (E) using Newton-Raphson method
        E = M  # Initial guess for E
        tolerance = 1e-10  # Set a tolerance for the solution
        while True:
            delta_E = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
            E -= delta_E
            if abs(delta_E) < tolerance:
                break
        
        # Calculate True Anomaly (ν) from Eccentric Anomaly (E)
        ν = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
        
        # Ensure ν is in the range [0, 2π]
        ν = np.mod(ν, 2 * np.pi)
        
        return np.degrees(ν)  # Convert ν to degrees if needed


    def find_velocity_from_elements(self):
        # Constants
        mu = G * M  # Gravitational parameter
        
        # True anomaly (ν) - assuming you have a method for this
        ν = self.calculate_true_anomaly()
        
        # Perifocal velocity components
        v_r = np.sqrt(mu / self.elements.a) * self.elements.eccentricity * np.sin(ν) / (1 + self.elements.eccentricity * np.cos(ν))
        v_t = np.sqrt(mu / self.elements.a) * (1 + self.elements.eccentricity * np.cos(ν)) / np.sqrt(1 - self.elements.eccentricity**2)
        
        # Rotate to geocentric equatorial frame
        Ω = np.radians(self.elements.raan)
        i = np.radians(self.elements.inclination)
        ω = np.radians(self.elements.arg_perigee)
        
        # Rotation matrix from perifocal to geocentric equatorial frame
        R = np.array([
            [np.cos(Ω) * np.cos(ω) - np.sin(Ω) * np.sin(ω) * np.cos(i), -np.cos(Ω) * np.sin(ω) - np.sin(Ω) * np.cos(ω) * np.cos(i), np.sin(Ω) * np.sin(i)],
            [np.sin(Ω) * np.cos(ω) + np.cos(Ω) * np.sin(ω) * np.cos(i), -np.sin(Ω) * np.sin(ω) + np.cos(Ω) * np.cos(ω) * np.cos(i), -np.cos(Ω) * np.sin(i)],
            [np.sin(ω) * np.sin(i), np.cos(ω) * np.sin(i), np.cos(i)]
        ])
        
        # Perifocal velocity vector
        v_perifocal = np.array([v_r * np.cos(ν) - v_t * np.sin(ν), v_r * np.sin(ν) + v_t * np.cos(ν), 0])
        
        # Velocity in geocentric equatorial frame
        velocity = np.dot(R, v_perifocal)
        vx, vy, vz = velocity
        
        return [vx, vy, vz]


    def find_elements_from_velocity(self):
        mu = G * M  # Standard gravitational parameter
        r_vec = np.array(self.position)  # Position vector
        v_vec = np.array(self.velocity)  # Velocity vector
        r = np.linalg.norm(r_vec)
        v = np.linalg.norm(v_vec)

        # Calculate specific angular momentum
        h_vec = np.cross(r_vec, v_vec)
        h = np.linalg.norm(h_vec)

        # Calculate eccentricity vector
        e_vec = np.cross(v_vec, h_vec) / mu - r_vec / r
        e = np.linalg.norm(e_vec)

        # Calculate energy
        epsilon = v**2 / 2 - mu / r

        # Semi-major axis
        a = -mu / (2 * epsilon)

        # Inclination
        i = np.arccos(h_vec[2] / h)

        # Node vector
        n_vec = np.cross([0, 0, 1], h_vec)
        n = np.linalg.norm(n_vec)

        # Right Ascension of the Ascending Node (RAAN)
        Omega = 0
        if n != 0:
            Omega = np.arccos(n_vec[0] / n)
            if n_vec[1] < 0:
                Omega = 2 * np.pi - Omega
        

        # Argument of Perigee
        omega = 0
        if n*e != 0:
            omega = np.arccos(np.dot(n_vec, e_vec) / (n * e))
            if e_vec[2] < 0:
                omega = 2 * np.pi - omega
        
        # True Anomaly
        nu = np.arccos(np.dot(e_vec, r_vec) / (e * r))
        if np.dot(r_vec, v_vec) < 0:
            nu = 2 * np.pi - nu

        # Convert radians to degrees for angles
        i_deg = np.degrees(i)
        Omega_deg = np.degrees(Omega)
        omega_deg = np.degrees(omega)
        nu_deg = np.degrees(nu)

        return orbital_elements(i_deg, Omega_deg, e, omega_deg, nu_deg, a)

    

    def find_position_from_elements(self):
        # Unpack orbital elements for readability
        i = np.radians(self.elements.inclination)
        Ω = np.radians(self.elements.raan)
        e = self.elements.eccentricity
        ω = np.radians(self.elements.arg_perigee)
        M = np.radians(self.elements.mean_anomaly)
        a = self.elements.a
        
        # Step 1: Solve Kepler's Equation for Eccentric Anomaly (E)
        E = solve_keplers_equation(M, e)
        
        # Step 2: Calculate True Anomaly (ν)
        ν = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
        
        # Step 3: Distance (r)
        r = a * (1 - e * np.cos(E))
        
        # Step 4: Position in the orbital plane (perifocal coordinates)
        x_p = r * np.cos(ν)
        y_p = r * np.sin(ν)
        
        # Step 5: Transform to the geocentric equatorial frame
        # Rotation matrix for argument of perigee (ω), inclination (i), and RAAN (Ω)
        R = np.array([
            [np.cos(Ω) * np.cos(ω) - np.sin(Ω) * np.sin(ω) * np.cos(i), -np.cos(Ω) * np.sin(ω) - np.sin(Ω) * np.cos(ω) * np.cos(i), np.sin(Ω) * np.sin(i)],
            [np.sin(Ω) * np.cos(ω) + np.cos(Ω) * np.sin(ω) * np.cos(i), -np.sin(Ω) * np.sin(ω) + np.cos(Ω) * np.cos(ω) * np.cos(i), -np.cos(Ω) * np.sin(i)],
            [np.sin(ω) * np.sin(i), np.cos(ω) * np.sin(i), np.cos(i)]
        ])
        
        position = np.dot(R, np.array([x_p, y_p, 0])) # Position in geocentric equatorial coordinates
        return position.tolist()
    

    def find_grav_acceleration(self):
        r = np.linalg.norm(self.position)
        r_hat = np.array(self.position) / r
        z = self.position[2]  # z-component of the position vector

        # Constants
        R_Earth = 1  # Mean radius of Earth, in meters
        

        # Gravitational acceleration without J2
        a_magnitude = G * M / (r**2)
        a_gravity = -a_magnitude * r_hat

        # J2 perturbations
        a_r_J2 = -1.5 * J2 * (G * M * R_Earth**2) / (r**4) * (1 - 5 * (z/r)**2)
        a_n_J2 = -1.5 * J2 * (G * M * R_Earth**2) / (r**4) * (3 * (z/r)**2)

        # Adding J2 perturbations to the radial component
        #a_J2 = a_r_J2 * r_hat + a_n_J2 * np.array([0, 0, 1])

        # Calculate the angular momentum vector
        h_vec = np.cross(self.position, self.velocity)
        h_hat = h_vec / np.linalg.norm(h_vec)

        # Normal component of J2 acceleration should be in the direction of h_hat
        a_n_J2 = -1.5 * J2 * (G * M * R_Earth**2) / (r**4) * (3 * (z/r)**2) * h_hat

        # Now combine the radial and normal components
        a_J2 = a_r_J2 * r_hat + a_n_J2






        # Total gravitational acceleration
        a_total = a_gravity + a_J2

        return a_total
    


    





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
        return f"arg_p = {self.arg_perigee}"
        

def solve_keplers_equation(M, e, tolerance=1e-6):
        E = M  # Initial guess
        delta = 1  # Initialize delta to enter the loop
        while abs(delta) > tolerance:
            delta = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
            E -= delta
        return E