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

# query JPL horizons for orbital elements and position vectors
obj = Horizons(id='Ceres', epochs={'start':'2010-01-01', 'stop':'2015-03-01', 'step':'100d'})
print(obj)

if os.path.isfile("df_el.csv") and os.path.isfile("df_vec.csv"):
    print("load files")
    df_el = pd.read_csv("df_el.csv", index_col = 0)
    df_vec = pd.read_csv("df_vec.csv", index_col = 0)

else:
    print("query JPL")
    el = obj.elements()
    vec = obj.vectors()
    df_el = el.to_pandas()
    df_vec = vec.to_pandas()
    df_el.to_csv("df_el.csv")
    df_vec.to_csv("df_vec.csv")

print(df_el)
print(df_vec)
print(list(df_el))

# Load the first set of orbital elements
Ceres = BodyOrb()
i=0
orb_dict = {"a":df_el.iloc[i]["a"],
"e":df_el.iloc[i]["e"],
"inc":np.radians(df_el.iloc[i]["incl"]),
"peri":np.radians(df_el.iloc[i]["w"]),
"node":np.radians(df_el.iloc[i]["Omega"]),
"M":np.radians(df_el.iloc[i]["M"])}
print(orb_dict)
Ceres.load_dict(orb_dict)
Ceres.calc_values()
Ceres.calc_orb_vectors()
Ceres.print_orb()

# obtain the x, y, z position for plotting orbit
df_pos = Ceres.planet_orbit()

# exit()

# make 3d orbit plot to compare the calculated orbit to the JPL positions vectors
fig = plt.figure()
ax = plt.axes(projection='3d')

ax.scatter3D(0,0,0,marker="+",c="k")
ax.scatter3D(df_vec["x"],df_vec["y"],df_vec["z"])
ax.set_box_aspect((np.ptp(df_vec["x"]), np.ptp(df_vec["y"]), np.ptp(df_vec["z"])))  # aspect ratio is 1:1:1 in data space

ax.plot3D(df_pos["x"],df_pos["y"],df_pos["z"])

plt.show()
