import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ensure tests can import the package
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from nice_orbs.orb_class import BodyOrb
import nice_orbs.orb_funcs as orb_funcs


print("test")

# make a blank body
Earth = BodyOrb()
print(Earth.a)
Earth.print_pos()

# set a value manually
Earth.a = 10
Earth.print_orb()

# load values from a dictionary
orb_dict = {"a":1.0, "e":1e-2, "inc":np.radians(5), "peri":np.radians(10), "node":np.radians(15), "f":np.radians(20)}
print(orb_dict)
Earth.load_dict(orb_dict)
Earth.print_orb()

print(Earth.node)
print(type(Earth.node))
print(np.cos(Earth.node))

# Calculate the orbit unit vectors
# Earth.ep = orb_pos.calc_ep_vec(Earth.node,Earth.peri,Earth.inc)
# Earth.eQ = orb_pos.calc_eQ_vec(Earth.node,Earth.peri,Earth.inc)
Earth.calc_orb_vectors()
print(Earth.ep)
print(Earth.eQ)

# calculate required values from supplied orbital elements
Earth.calc_values()

# find pos and vel as a function of f
# f_true = np.linspace(0.0,2.0*np.pi)
# f_true = np.reshape(f_true,(len(f_true),1))
f_true = 0.5
df_rv = Earth.pos_vel_from_orbit(f_true)
print(df_rv)

df_pos = Earth.planet_orbit()
print(df_pos)

fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])

ax1.scatter(df_pos["x"],df_pos["y"],c=df_pos["f"])
ax1.plot(df_pos["x"],df_pos["y"])

ax1.scatter(df_rv["x"],df_rv["y"],c="k",marker="x")

ax1.set_aspect("equal")

plt.show()
