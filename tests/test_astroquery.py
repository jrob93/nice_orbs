"""
Test the xyz positions calculated by nice_orbs against the 'ground truth' positions from JPL horizons
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits import mplot3d
from astroquery.jplhorizons import Horizons
from astropy.table import Table

# ensure tests can import the package
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from nice_orbs.orb_class import BodyOrb
import nice_orbs.orb_funcs as orb_funcs

# set up 3d orbit plot to compare the calculated orbit to the JPL positions vectors
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(0,0,0,marker="+",c="k") # plot origin (sun)

# Go through different objects (3 = Earth, 90004507 = C/2020 H5)
body_list = [3,"Ceres","Hilda","Vesta","Don Quixote"]#,90004507]
body_type = [None,"name","name","name","name"]#,None]
for Bodyname, Bodytype in zip(body_list, body_type):

    el_name = "test_data/df_el_{}.csv".format(Bodyname)
    vec_name = "test_data/df_vec_{}.csv".format(Bodyname)

    # query JPL horizons for orbital elements and position vectors
    obj = Horizons(id=Bodyname, epochs={'start':'2010-01-01', 'stop':'2015-03-01', 'step':'100d'}, #id_type = Bodytype
    )
    print(obj)

    # the JPL queries will be saved and loaded to save re-run time
    if os.path.isfile(el_name) and os.path.isfile(vec_name):
        print("load files")
        df_el = pd.read_csv(el_name, index_col = 0)
        df_vec = pd.read_csv(vec_name, index_col = 0)

    else:
        print("query JPL")
        el = obj.elements()
        vec = obj.vectors()
        df_el = el.to_pandas()
        df_vec = vec.to_pandas()
        df_el.to_csv(el_name)
        df_vec.to_csv(vec_name)

    print(df_el)
    print(df_vec)
    print(list(df_el))

    # Load the first set of orbital elements
    Body = BodyOrb()
    i=0
    orb_dict = {"a":df_el.iloc[i]["a"],
    "e":df_el.iloc[i]["e"],
    "inc":np.radians(df_el.iloc[i]["incl"]),
    "peri":np.radians(df_el.iloc[i]["w"]),
    "node":np.radians(df_el.iloc[i]["Omega"]),
    "M":np.radians(df_el.iloc[i]["M"])}
    print(orb_dict)
    Body.load_dict(orb_dict)
    Body.calc_values()
    Body.calc_orb_vectors()
    Body.print_orb()

    # obtain the x, y, z position for plotting orbit using nice_orbs
    df_pos = Body.planet_orbit()

    # exit()

    ax.scatter3D(df_vec["x"],df_vec["y"],df_vec["z"],label = Bodyname)
    ax.plot3D(df_pos["x"],df_pos["y"],df_pos["z"])

ax.set_box_aspect((np.ptp(df_vec["x"]), np.ptp(df_vec["y"]), np.ptp(df_vec["z"])))  # aspect ratio is 1:1:1 in data space
# ax.set_box_aspect((10,10,10))

ax.legend()
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()
