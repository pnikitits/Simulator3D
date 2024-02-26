import numpy as np
from Constants import *
from Extra import clean_3vec
from astropy import units as u
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from Orbital_Elements import *


class Satellite:
    def __init__(self , elements):
        self.elements = elements
        self.position , self.velocity = self.find_pos_vel()
        
        self.position = clean_3vec(self.position)
        self.velocity = clean_3vec(self.velocity)

        # print(f"\nvel = {self.velocity}")
        # print(f"pos = {self.position}\n")
        

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
        return float(np.sqrt((G*M) / (self.elements.a**3))) # [rad / s]


