import numpy as np
from Constants import *



def normalize_vector(vector):
    norm = np.linalg.norm(vector)
    if norm == 0:
        return vector # Or handle this case as needed, maybe return None or raise an error
    return vector / norm


def angle_between_vectors(pos1, pos2 , vel1 , deg=True):

    normal_axis = np.cross(pos1 , vel1)
    dot_product = np.dot(pos1, pos2)
    magnitude_v1 = np.linalg.norm(pos1)
    magnitude_v2 = np.linalg.norm(pos2)
    cos_angle = dot_product / (magnitude_v1 * magnitude_v2)
    angle_radians = np.arccos(np.clip(cos_angle, -1.0, 1.0))
    cross_product = np.cross(pos1, pos2)
    direction = np.dot(cross_product, normal_axis)
    
    if direction < 0:
        angle_radians = 2*np.pi - angle_radians
    if deg:
        return np.degrees(angle_radians)
    else:
        return angle_radians
    

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


def clean_3vec(x):
    return [float(x[0]) , float(x[1]) , float(x[2])]

