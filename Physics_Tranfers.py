import numpy as np
from Constants import *
from Satellite import *
from Extra import *


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




def algorithm_45(otv:Satellite , target:Satellite , prints=False , debug_msg=False):
    """
    Modified Algorithm 45 (page 363), diagram (page 361)

    Structure:
        1. Compute the angular velocities [rad/s]
        2. Compute Hohmann transfer time
        3. Compute lead, phase and initial phase angles:
            - we need to trigger 1st Hohmann boost when phase angle == initial phase angle
            - if we are past the phase angle, += 360 to the angle to go
        4. Compute the Phasing time (to wait before dv1)

    In:
        otv       : Satellite
        target    : Satellite
        prints    : Bool (for the HUD)
        debug_msg : Bool

    Return:
        T_wait | lead_angle , phase_angle , initial_phase_angle , T_wait
    """

    print(f"\n----- algorithm_45 start -----\n") if debug_msg else None

    target_a = target.elements.a # [m]
    otv_a = otv.elements.a       # [m]

    mu = G*M
    omega_tgt = np.sqrt(mu / target_a**3)
    omega_otv = np.sqrt(mu / otv_a**3)
    print(f"omega_tgt = {omega_tgt}\nomega_otv = {omega_otv}") if debug_msg else None

    T_transfer = hohmann_time(otv_a , target_a)
    print(f"T_transfer = {T_transfer}") if debug_msg else None

    lead_angle = omega_tgt * T_transfer
    print(f"lead_angle = {np.degrees(lead_angle)}") if debug_msg else None
    
    phase_angle = abs(lead_angle - np.pi) # This must be positive
    print(f"phase_angle = {np.degrees(phase_angle)}") if debug_msg else None

    print(otv.position , target.position , otv.velocity) if debug_msg else None
    initial_phase_angle = angle_between_vectors(otv.position , target.position , otv.velocity , deg=False)
    print(f"initial_phase_angle = {(initial_phase_angle)}") if debug_msg else None

    angle_to_go = initial_phase_angle - phase_angle
    if angle_to_go < 0:
        angle_to_go += 2*np.pi
    print(f"angle_to_go = {np.degrees(angle_to_go)}") if debug_msg else None
 
    T_wait = angle_to_go / (omega_otv - omega_tgt)
    print(f"T_wait = {T_wait}") if debug_msg else None

    print(f"\n----- algorithm_45 end -----\n") if debug_msg else None
    if prints:
        return lead_angle , phase_angle , initial_phase_angle , T_wait
    return T_wait



def simple_phase(object:Satellite , target_anomaly):
    """
    Calculates phasing time to reach a static point on the orbit
    """
    d_ma = np.radians(target_anomaly - object.elements.mean_anomaly_360)

    if d_ma < 0:
        d_ma += 2*np.pi
        
    phase_t =  d_ma / object.find_angular_velocity()

    return phase_t / DT # time in seconds