import numpy as np
import pandas as pd
from . import orb_funcs

# Define the orbiting body class in this package

class BodyOrb:

    def __init__(self):

        # time/epoch
        self.t = None

        # cartesian x, y, z position
        self.x = None
        self.y = None
        self.z = None
        # cartesian x, y, z velocity
        self.vx = None
        self.vy = None
        self.vz = None

        # radial distance
        self.r = None

        # Keplerian elements
        self.a = None # semimajor axis
        self.e = None # eccentrcity
        self.inc = None # inclination
        self.peri = None # argument of pericentre
        self.node = None # longitude of ascending node
        self.f = None # true anomaly
        self.M = None # mean anomaly
        self.E = None # eccentric anomaly
        self.n = None # mean motion
        self.eta = None
        self.long_peri = None # longitude of pericentre = peri + node
        self.mean_long = None # mean longitude = peri + node + M
        self.peri_time = None # time of pericentre passage

        # Orbit vectors
        self.ep = None
        self.eQ = None

        self.G = 1.0 # Gravitational constant
        self.M = 1.0 # central body mass
        self.m = 0.0 # body mass if required
        self.mu = self.G*(self.M+self.m) # Gravitational parameter where M is central body mass


    def print_pos(self):
        """ print the positional elements """
        print("x={},y={},z={}".format(self.x,self.y,self.z))

    def print_orb(self):
        """ print the orbital elements """
        print("a={},e={},inc={},peri={},node={},f={}".format(
            self.a,self.e,self.inc,self.peri,self.node,self.f))

    def load_dict(self, x):
        """ load body parameters from a dictionary

        Parameters
        ----------
        x
            A dict of BodyOrb attributes, e.g. x = {"a":1.0, "e":1e-2, "inc":np.radians(5), "peri":np.radians(10), "node":np.radians(15), "f":np.radians(20)}
        """
        for y in x:
            print(y)
            exec("self.{} = {}".format(y, x[y])) # https://stackoverflow.com/questions/8307612/how-to-create-class-variable-dynamically-in-python

    def calc_orb_vectors(self):
        """ Find the unit vectors describing the orbit orientation
        Parameters
        ----------
        node
            longitude of ascending node
        peri
            argument of pericentre
        inc
            inclination
        """

        self.ep = orb_funcs.calc_ep_vec(self.node,self.peri,self.inc)
        self.eQ = orb_funcs.calc_eQ_vec(self.node,self.peri,self.inc)

    def calc_values(self):
        # !!! write a function to calculate missing orbital elements if only some are provided?

        self.n=np.sqrt(self.mu/(self.a**3))
        self.eta=np.sqrt(1.0-(self.e**2))

    def planet_orbit(self,n = 100):
        '''
        Function to find the xyz points which describe an orbit relative to the reference point (typically heliocentric)

        Parameters
        ----------
        self
            the BodyOrb class
        n
            Number of points used to calculate orbits
        '''

        #specify f_true from 0 to 2pi radians: i.e. number of points on orbit, THE TRUE ANOMALY
        # by going from 0 to exactly 2pi the first and last position will be the same so that a line plot will be a closed loop
        f_true = np.linspace(0.0,2.0*np.pi, n)
        f_true = np.reshape(f_true,(n,1))

        # find the radial distance of the orbit at all f
        r=orb_funcs.r_elliptical(self.a,self.e,f_true)

        # calculate the r(x,y,z) position array
        pos=r*((np.cos(f_true)*self.ep)+(np.sin(f_true)*self.eQ)) # PSS eq 11.36a

        # return a dataframe of true anomaly and x, y, z position
        df_pos = pd.DataFrame(data = np.hstack((f_true,pos)), columns = ["f","x","y","z"])
        return df_pos

    def pos_vel_from_orbit(self,f_true):
        '''function to find xyz and vxvyvz from the a,e,I,OMEGA,omega,f_true orbit data, where mu=G(M+m)


        Parameters
        ----------
        self
        f_true
            value (or array) of true anomaly. If array must be correct shape - (N,1)
        '''

        # calculate position and velocity arrays and combine with f
        r = orb_funcs.r_elliptical(self.a,self.e,f_true)
        pos=r*((np.cos(f_true)*self.ep)+(np.sin(f_true)*self.eQ))#PSS eq 11.36a
        vel=(self.n*self.a/self.eta)*((-np.sin(f_true)*self.ep)+((self.e+np.cos(f_true))*self.eQ))
        data = np.hstack((np.hstack((f_true,pos)),vel))

        # if f_true is just a single value we must reshape
        if not isinstance(f_true, np.ndarray):
            data = np.reshape(data,(1,7))

        # return a dataframe with position and velocity
        df_pos_vel = pd.DataFrame(data, columns = ["f","x","y","z","vx","vy","vz"])
        return df_pos_vel

    # pos and vel as a function of other parameters, e.g. time, mean anomaly etc

    # Cartesian position and velocity to orbit
