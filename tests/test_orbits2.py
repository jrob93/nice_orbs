import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ensure tests can import the package
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from nice_orbs.orb_class import BodyOrb
import nice_orbs.orb_funcs as orb_funcs

# get eccentric anomaly from a mean anomaly
e = 0.7
M = np.radians(45)
E = orb_funcs.E_from_M(M,e)
print(E)
