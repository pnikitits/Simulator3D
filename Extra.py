import numpy as np
G = 0.01 #6.67e-11 # 
M = 2000 #5.972e24 # 



def move_object(object , xyz):
    x = object.getX()
    y = object.getY()
    z = object.getZ()

    object.setX(x + xyz[0])
    object.setY(y + xyz[1])
    object.setZ(z + xyz[2])
        

def rotate_object(object , xyz):
    x = object.getH()
    y = object.getP()
    z = object.getR()

    object.setH(x + xyz[0])
    object.setP(y + xyz[1])
    object.setR(z + xyz[2])



def polar_to_cartesian(angle, radius):
    # Convert angle to radians
    angle_rad = np.radians(angle)
    # Calculate Cartesian coordinates
    x = radius * np.cos(angle_rad)
    y = radius * np.sin(angle_rad)
    return [x, y]


def find_circular_orbit_v(central_obj , rad):
    return np.sqrt(G * central_obj.mass / rad)

def normalize_vector(vector):
    return vector / np.linalg.norm(vector)

def normalize_vector(vector):
    norm = np.linalg.norm(vector)
    if norm == 0:
        return vector  # Or handle this case as needed, maybe return None or raise an error
    return vector / norm

def calculate_distance(v1 , v2):
    return np.linalg.norm(np.array(v2) - np.array(v1))

def calculate_vector(v1 , v2):
    return np.array(v2) - np.array(v1)



def calc_semi_major_axis(v , r , c_obj , G=G):
    """
    In:
        v : current velocity
        r : current orbit radius (not altitude)
        c_obj : centre object
    Out:
        a : semi-major axis 
    """
    M = c_obj.mass
    v = np.linalg.norm(v)

    p1 = -G*M
    p2 = 2*( (v**2)/2 - G*M/r )

    a = p1/p2
    return a


def calc_eccentricity(v , r , c_obj , G=G):
    """
    In:
        v : current velocity
        r : current orbit radius (not altitude)
        c_obj : centre object
    Out:
        e : eccentricity
    """
    M = c_obj.mass
    v = np.linalg.norm(v)

    p1 = 2*( (v**2)/2 - G*M/r ) * (r*v)**2
    p2 = (G*M)**2

    e = np.sqrt(1 + p1/p2)
    return e


def calc_semi_minor_axis(a , e):
    """
    In:
        a : semi-major axis
        e : eccentricity
    Out:
        b : semi-minor axis
    """
    if e > 1:
        e = 1

    b = a * np.sqrt(1 - e**2)
    return b



def average_over_chunks_with_none(values , stride=5):
    out_ = []

    for i, value in enumerate(values):
        if i % stride == 0:
            out_.append(np.mean(values[i:i+stride]))
        else:
            out_.append(0)
    print(len(out_)) 
    out_ = list(out_)
    print(out_)

    
    for i in range(0 , len(out_) , stride):
        sub_rg = np.linspace(out_[i] , out_[i+stride-1] , stride)
        for j in range(0 , stride-1):
            out_[i+j] = sub_rg[j]

    print(out_)


    return out_
    



