import numpy as np
from Extra import G , M
import math
import random
from datetime import datetime, timedelta


class Debris:
    def __init__(self , norad , inclination , raan , eccentricity ,
                 arg_perigee , mean_anomaly , a , rcs):
        
        # ID
        self.norad = norad

        # Angle between orbital plane and Earth equator plane
        self.inclination = inclination

        # Right ascension of ascending node (orientation of orbital plane)
        self.raan = raan

        # 0 for circular orbit, 0<e<1 for ellipic orbit
        self.eccentricity = eccentricity

        # Location of closest approach wrt ascending node
        self.arg_perigee = arg_perigee

        # Angular position in orbit (from perigee)
        self.mean_anomaly = mean_anomaly

        # Semi-major axis
        self.a = a

        # Radar cross section
        self.rcs = rcs

        time = calculate_time()

        # self.initial_pos = self.calculate_position(central_mass=M , grav_const=G , current_time=time)
        # self.initial_vel = self.calculate_velocity(central_mass=M , grav_const=G , current_time=time)




    # def calculate_position(self, central_mass, grav_const, current_time):
    #     """
    #     Calculate the position of the debris at a given time

    #     Parameters:
    #     - central_mass: Mass of the central object
    #     - grav_const: Gravitational constant
    #     - current_time: Current time in seconds

    #     Returns:
    #     x, y, z coordinates of the debris
    #     """
    #     # Calculate mean motion
    #     n = math.sqrt(grav_const * central_mass / math.pow(self.a, 3))

    #     # Calculate mean anomaly at current time
    #     mean_anomaly_at_time = self.mean_anomaly + n * current_time

    #     # Calculate eccentric anomaly using iterative method (e.g., Newton's method)
    #     eccentric_anomaly = self.solve_keplers_equation(mean_anomaly_at_time)

    #     # Calculate true anomaly
    #     true_anomaly = 2 * math.atan(math.sqrt((1 + self.eccentricity) / (1 - self.eccentricity)) * math.tan(eccentric_anomaly / 2))

    #     # Calculate distance from central object
    #     r = self.a * (1 - self.eccentricity ** 2) / (1 + self.eccentricity * math.cos(true_anomaly))

    #     # Calculate position coordinates
    #     x = r * (math.cos(self.raan) * math.cos(true_anomaly + self.arg_perigee) - math.sin(self.raan) * math.sin(true_anomaly + self.arg_perigee) * math.cos(self.inclination))
    #     y = r * (math.sin(self.raan) * math.cos(true_anomaly + self.arg_perigee) + math.cos(self.raan) * math.sin(true_anomaly + self.arg_perigee) * math.cos(self.inclination))
    #     z = r * (math.sin(true_anomaly + self.arg_perigee) * math.sin(self.inclination))

    #     return np.array([x, y, z])

    # def calculate_velocity(self, central_mass, grav_const, current_time):
    #     """
    #     Calculate the velocity of the debris at a given time

    #     Parameters:
    #     - central_mass: Mass of the central object
    #     - grav_const: Gravitational constant
    #     - current_time: Current time in seconds

    #     Returns:
    #     x, y, z components of the velocity
    #     """
    #     # mean motion
    #     n = math.sqrt(grav_const * central_mass / math.pow(self.a, 3))

    #     # mean anomaly at current time
    #     mean_anomaly_at_time = self.mean_anomaly + n * current_time

    #     # eccentric anomaly
    #     eccentric_anomaly = self.solve_keplers_equation(mean_anomaly_at_time)

    #     # true anomaly
    #     true_anomaly = 2 * math.atan(math.sqrt((1 + self.eccentricity) / (1 - self.eccentricity)) * math.tan(eccentric_anomaly / 2))

    #     # specific angular momentum
    #     h = math.sqrt(grav_const * central_mass * self.a * (1 - self.eccentricity ** 2))

    #     # velocity components
    #     vx = (h / self.a) * (-math.sin(true_anomaly))
    #     vy = (h / self.a) * (self.eccentricity + math.cos(true_anomaly))
    #     vz = (h / self.a) * math.sin(true_anomaly) * math.sin(self.inclination)

    #     return np.array([vx, vy, vz])


    # def solve_keplers_equation(self, mean_anomaly):
    #     """
    #     Solve Kepler's equation for eccentric anomaly using Newton's method.

    #     Parameters:
    #     - mean_anomaly: Mean anomaly in radians

    #     Returns:
    #     Eccentric anomaly in radians
    #     """
    #     # Initial guess for eccentric anomaly
    #     E = mean_anomaly

    #     # Iterative solution using Newton's method
    #     max_iterations = 100
    #     tolerance = 1e-10

    #     for _ in range(max_iterations):
    #         E_next = E - (E - self.eccentricity * math.sin(E) - mean_anomaly) / (1 - self.eccentricity * math.cos(E))

    #         # Check for convergence
    #         if abs(E_next - E) < tolerance:
    #             break

    #         E = E_next

    #     return E
    

def make_random_debris(n):
    output = []
    for i in range(n):
        # Make some random values
        inclination = random.uniform(0 , np.pi)
        raan = random.uniform(0 , 2*np.pi)
        eccentricity = 0.0#random.uniform(0 , 1)
        arg_perigee = random.uniform(0 , 2*np.pi)
        mean_anomaly = random.uniform(0 , 2*np.pi)
        a = random.uniform(1.2 , 2.0)

        debris = Debris(norad=f"Debris_{i}",
                        inclination=inclination,
                        raan=raan,
                        eccentricity=eccentricity,
                        arg_perigee=arg_perigee,
                        mean_anomaly=mean_anomaly,
                        a=a,
                        rcs=None)
        output.append(debris)

    return output


# def calculate_time():
    

#     # Define the specific date
#     specific_date = datetime(2017, 5, 6, 0, 0, 0)  # May 6, 2017, 00:00:00

#     # Get the current date and time
#     current_date = datetime.utcnow()

#     # Calculate the time difference
#     time_difference = current_date - specific_date

#     # Convert the time difference to seconds
#     time_difference_seconds = time_difference.total_seconds()

#     return time_difference_seconds

