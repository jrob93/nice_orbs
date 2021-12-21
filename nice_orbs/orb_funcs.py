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
