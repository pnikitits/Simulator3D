import numpy as np
from Constants import *
from Satellite import *


def hohmann_dv(r1 , r2):
    # Indepent of moving up or down
    mu = G*M

    dv1 = np.sqrt(mu/r1)*(np.sqrt(2*r2/(r1+r2))-1)
    dv2 = np.sqrt(mu/r2)*(np.sqrt(2*r1/(r1+r2))-1)

    if r1 < r2:
        return abs(dv1) , abs(dv2)
    elif r1 > r2:
        return -abs(dv1) , -abs(dv2)



def hohmann_time(r1 , r2):
    # Indepent of moving up or down
    mu = G*M
    th = np.pi * np.sqrt( ((r1+r2)**3) / (8*mu) )
    return th



def inc_dv(v1 , inc1 , inc2):
    inc1 = np.deg2rad(inc1)
    inc2 = np.deg2rad(inc2)

    return 2*v1*np.sin((inc2-inc1)/2)


def raan_dv(v1 , raan1 , raan2):
    raan1 = np.deg2rad(raan1)
    raan2 = np.deg2rad(raan2)

    return 2*v1*np.sin((raan2-raan1)/2)


def inc_and_raan_dv(v1 , raan1 , raan2 , inc1 , inc2):
    raan1 = np.deg2rad(raan1)
    raan2 = np.deg2rad(raan2)
    inc1 = np.deg2rad(inc1)
    inc2 = np.deg2rad(inc2)

    d_raan = raan2 - raan1
    theta = np.arccos( np.cos(d_raan)*np.sin(inc1)*np.sin(inc2) + np.cos(inc1)*np.cos(inc2) )
    return 2*v1*np.sin(theta/2)




def delta_u(r1 , r2):
    return np.pi * (1 - np.sqrt( (r1+r2) / (2*r2) ))


def phase_time(otv:Satellite , target:Satellite):
    
    r1 = np.linalg.norm(otv.position)
    r2 = np.linalg.norm(target.position)
    ang_1 = np.deg2rad(otv.elements.mean_anomaly)
    ang_2 = np.deg2rad(target.elements.mean_anomaly)

    # Check orbits
    if r1 == r2:
        return 0
    
    # delta_u
    du = delta_u(r1 , r2)
    if r1 > r2:
        reverse_du = delta_u(r2 , r1)
        du = reverse_du

    # angle diff
    angle_diff = ang_2 - ang_1
    if r1 < r2 and angle_diff < du:
        angle_diff += 2*np.pi
    elif r1 > r2 and angle_diff > du:
        angle_diff -= 2*np.pi
    
    # angle velocity
    ang_vel_1 = np.sqrt(G*M / (r1**3))
    ang_vel_2 = np.sqrt(G*M / (r2**3))
    
    # phasing time
    dt = (du - angle_diff) / (ang_vel_2 - ang_vel_1)
    
    return dt



def algorithm_45(otv:Satellite , target:Satellite):

    mu = G*M
    omega_tgt = np.sqrt(mu / target.elements.a**3)
    omega_otv = np.sqrt(mu / otv.elements.a**3)

    T_transfer = hohmann_time(otv.elements.a , target.elements.a)
    lead_angle = omega_tgt * T_transfer
    
    phase_angle = lead_angle - np.pi # This must be positive
    initial_phase_angle = np.radians(target.elements.mean_anomaly - otv.elements.mean_anomaly)

    T_wait = (phase_angle - initial_phase_angle + 2*np.pi*1) / (omega_otv - omega_tgt)
    return T_wait



if __name__ == "__main__":
    sat_1 = Satellite(orbital_elements(inclination=0,
                                       raan=0,
                                       eccentricity=0,
                                       arg_perigee=0,
                                       mean_anomaly=0,
                                       a=2))
    
    sat_2 = Satellite(orbital_elements(inclination=0,
                                       raan=0,
                                       eccentricity=0,
                                       arg_perigee=0,
                                       mean_anomaly=0,
                                       a=4))
    
    print(algorithm_45(otv=sat_1 , target=sat_2))