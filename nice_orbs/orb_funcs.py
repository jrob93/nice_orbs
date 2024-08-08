import numpy as np

# These functions to calculate orbital values and vectors can be used separately, or on the BodyOrb class

def r_elliptical(a,e,f):
    '''
    Function to find distance, r, at true anomaly, f, around a keplerian orbit

    Parameters
    ----------
    a
        semi major axis (m)
    e
        orbital eccentricity
    f
        true anomaly (rad)
    '''
    r = (a*(1.0-e**2))/(1.0+e*np.cos(f))
    return r

def calc_ep_vec(node,peri,inc):
    '''
    Function to find the normalised Lenz vector along the apsidal line, uses PSS eq11.39a

    Parameters
    ----------
    node
        longitude of ascending node
    peri
        argument of pericentre
    inc
        inclination
    '''
    epx=(np.cos(node)*np.cos(peri))-(np.cos(inc)*np.sin(node)*np.sin(peri)) # x component
    epy=(np.sin(node)*np.cos(peri))+(np.cos(inc)*np.cos(node)*np.sin(peri)) # y component
    epz=np.sin(inc)*np.sin(peri) # z component

    return  np.array([epx, epy, epz])

def calc_eQ_vec(node,peri,inc):
    '''
    Function to find the normalised e_Q vector which is perp to h and e_a and
    lies in the orbital plane, uses PSS eq11.39b

    Parameters
    ----------
    node
        longitude of ascending node
    peri
        argument of pericentre
    inc
        inclination
    '''
    eQx=-(np.cos(node)*np.sin(peri))-(np.cos(inc)*np.sin(node)*np.cos(peri)) #x
    eQy=-(np.sin(node)*np.sin(peri))+(np.cos(inc)*np.cos(node)*np.cos(peri)) #y
    eQz=np.sin(inc)*np.cos(peri) #z

    return np.array([eQx, eQy, eQz])
    
def M_from_E(E,e):
    """
    Find mean anomaly from the eccentric anomaly

    Parameters
    ----------
    E
        eccentric anomaly
    e
        eccentricity
    """
    M = E - (e*np.sin(E))
    return M

def M_from_t(t,n,M0,t0):
    """
    find mean anomaly at time t

    Parameters
    ----------
    t
        time to find M for
    n
        mean motion
    M0
        reference mean anomaly at time t0
    t0
        reference time t0

    Note that M = 0 at pericentre so if t0 = pericentre passage time then M0 = 0
    """

    M = M0 + (n*(t-t0))
    return M

def E_from_M(M,e,E_limit=1e-6):
    '''
    Solve Kepler's Equation, to find E from M, to within a certain limit.
    This is an iterative method and DIVERGES for e>0.662743
    See Murray and Dermott Solar System Dynamics p.35

    Parameters
    ----------
    M
        mean anomaly, right now only works for a single value of M, not an array
    e
        eccentricity
    E_limit
        Fractional cut off for the change in E between iterations
    '''

    if e>0.662743: # does this catch arrays too?
        # !!! throw a proper warning
        print("WARNING: high e, beware of divergence!")

    E_0=M
    deltaE=2.0*E_limit #set deltaE to some value, greater than the limit
    while deltaE>E_limit: # iterate until the change in E is less than the selected limit
        E_1=M+e*np.sin(E_0)
        deltaE=np.absolute(E_1-E_0)
        E_0=E_1
    return E_0

# !!! Add another function that is more accurate for high eccentrcity orbits?

def f_from_E(E,e):
    '''
    Function to find true anomaly f from eccentric anomaly E

    Parameters
    ----------
    E
        Eccentric anomaly
    e
        eccentricity
    '''

    f_true=2.0*np.arctan(np.sqrt((1.0+e)/(1.0-e))*np.tan(E/2.0))
    return f_true

def f_from_M(M,e):
    """
    Find true anomaly f from mean anomaly M

    Parameters
    ----------
    M
        mean anomaly
    e
        eccentricity
    """

    f_true=2.0*np.arctan(np.sqrt((1.0+e)/(1.0-e))*np.tan(E_from_M(M,e)/2.0))
    return f_true

def f_from_t(t,e,n,M0,t0):
    """
    find true anomaly f from time t

    Parameters
    ----------
    t
        time
    e
        eccentricity
    n
        mean motion
    M0
        reference mean anomaly
    t0
        reference time
    """
    M = M_from_t(t,n,M0,t0)
    f_true = f_from_M(M,e)
    return f_true

def rot_matrix(alpha = 0, beta = 0, gamma = 0, order = [0, 1, 2]):
    """Rotation about the x, y, z axes (right handed) in a given order. All angles are in radians.
    alpha:
        angle of right handed rotation around x axis using matrix Rx
    beta:
        angle of right handed rotation around y axis using matrix Ry
    gamma:
        angle of right handed rotation around z axis using matrix Rz
    order:
        order of operation for rotation about each axis, where Rx = 0, Ry = 1, Rz = 2.
        E.g. for rotating around x, then y, then z: R = Rz . Ry . Rx and order = [0,1,2]
        For rotating around z, then y, then x: R = Rx . Ry . Rz and order = [2,1,0]

    """
    
    Rx = np.array([[1, 0, 0],[0, np.cos(alpha), -np.sin(alpha)],[0, np.sin(alpha), np.cos(alpha)]])
    Ry = np.array([[np.cos(beta), 0, np.sin(beta)],[0, 1, 0],[-np.sin(beta), 0, np.cos(beta)]])
    Rz = np.array([[np.cos(gamma), -np.sin(gamma), 0],[np.sin(gamma), np.cos(gamma), 0],[0, 0, 1]])
    
    R_list = [Rx, Ry, Rz]
    
    R = np.dot(np.dot(R_list[order[2]], R_list[order[1]]), R_list[order[0]]) # the rotation order is important!
    
    return R

def mutual_ascending_node_f_true(bod1,bod2):
    """
    Return the f_true (in bod1 and bod2 elements) for the mutual ascending node of bod2 relative to bod1.
    N.B. to get the descending node, use node_mutual_f2 += np.pi
    
    bod1:
        The reference orbit
    bod2:
        Orbit that we are determining mutual ascending node for
        
    Returns
    
    node_mutual_f2,node_mutual_f2
    """
    
    # define the transforms
    trans_bod1 = rot_matrix(alpha = bod1.inc, gamma = bod1.node) # transform from cartesian to bod1 coplanar frame
    _trans_bod1 = np.linalg.inv(trans_bod1) # inverse: transform from bod1 frame to cartesian
    
    # define the required unit vectors and their transforms
    unit_x = np.array([1,0,0])
    unit_z = np.array([0,0,1])
    _unit_x = np.dot(trans_bod1,unit_x)
    _unit_z = np.dot(trans_bod1,unit_z)
    
    # transform the bod2 orbital vectors to the bod1 coplanar reference frame (?)
    bod2_ep = bod2.ep
    bod2_eQ = bod2.eQ
    _bod2_ep = np.dot(_trans_bod1, bod2_ep)
    _bod2_eQ = np.dot(_trans_bod1, bod2_eQ)
        
    # determine the f_true of the mutual ascending node for bod2
    # calculated when z component of position vector (in the bod1 reference frame) equals zero
    node_mutual_f2 = np.arctan2(-_bod2_ep[2],_bod2_eQ[2])
    # print("mutual node f (bod2) = {}".format(np.degrees(node_mutual_f2)))
    
    # Perform checks on the angle:
    # Ensure that this is the ascending node
    _node2 = bod2.r_vec(node_mutual_f2)
    _node_vel2 = bod2.v_vec(node_mutual_f2)
    # if np.dot(_node_vel2,_unit_z)>0:
    #     # print("ascending node")
    if np.dot(_node_vel2,_unit_z)<0:
        # print("descending node, convert to ascending")
        node_mutual_f2 += np.pi # get the other node
    # else:
    #     print("node ill defined?")
    #     print(np.dot(_node_vel2,_unit_z))
    # Also ensure node_mutual_f2 is [0,2pi]
    while node_mutual_f2 < 0: 
        # print("correct node_mutual_f2")
        node_mutual_f2 += (2.0*np.pi) 

    # determine the mutual nodal angle for bod1
    # this is the angle between the bod1 frame x vector and the bod2 nodal vector (which was determined in the bod1 ref frame)
    node_mutual1 = np.arctan2(np.linalg.norm(np.cross(_unit_x,_node2)),np.dot(_unit_x,_node2))
    # N.B. node_mutual1 is the smallest angle between the vectors, not the right hand rotation in f_true,
    # correct the angle as required so that the rotation is correct
    if (np.dot(_trans_bod1,np.cross(_unit_x,_node2))[2] < 0) and (node_mutual1<np.pi):
        # print("correct node_mutual1")
        node_mutual1 = (2.0*np.pi) - node_mutual1
    # convert this nodal angle to f_true in the bod1 elements by accounting for the argument of periapsis
    node_mutual_f1 = node_mutual1 - bod1.peri
    # print("mutual node f (bod1) = {}".format(np.degrees(node_mutual_f1)))

    return node_mutual_f1,node_mutual_f2