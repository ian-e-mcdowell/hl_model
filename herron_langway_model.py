# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 16:08:51 2022

@author: imcdo
"""

''' 
This program calculates the steady state Herron and Langway firn density 
profile given average temperature, accumulation rate, and initial surface snow
density.

INPUT VARIABLES:
Tavg: 10 m temperature ( or average site temperature) in celsius
bdot: accumulation rate in mwe/yr or (kg/m2/yr)
rhos: surface density in kg/m3
z: depth in meters

MODEL OUTPUT: 
rho: density (kg/m3) for all z-values.
zieq: ice equivalent depth for all z-values.
t: age for all z-values (only taking densification into account)

'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# Below here I just set up some matplotlib figure characteristics that will be applied for all figures created in this script.
# Font style and size
mpl.rcParams['font.family'] = 'Arial'         # Font
mpl.rcParams['font.size'] = 22                # General font size unless set below
mpl.rcParams['axes.labelsize'] = 24           # Axes labels font size
mpl.rcParams['figure.titlesize'] = 26         # Title font size
mpl.rcParams['figure.titleweight'] = 'bold'   # Bold title
mpl.rcParams['axes.labelweight'] = 'bold'     # Bold axes labels

# Axes and ticks parameters
mpl.rcParams['axes.linewidth'] = 1            # Width of axes border
mpl.rcParams['xtick.direction'] = 'inout'     # Make x ticks inside and outside of axis line
mpl.rcParams['ytick.direction'] = 'inout'     # Make y ticks inside and outside of axis line
mpl.rcParams['xtick.major.size'] = 12        # Set x tick length 
mpl.rcParams['ytick.major.size'] = 12         # Set y tick length
mpl.rcParams['xtick.major.width'] = 2         # Set x tick width 
mpl.rcParams['ytick.major.width'] = 2         # Set y tick width



# =============================================================================
# Define certain parameters used in the H&L model - don't change these
# =============================================================================
rhoi = 917     # density of ice (kg/m3)
rhoc = 550     # critical density where model equations change (kg/m3)
rhow = 1000    # density of water (kg/m3)
R = 8.314      # gas constant (J/(K mol))


# =============================================================================
# Input variables -- change these to your specific case
# =============================================================================
Tavg = -20 # set to mean annual air temp (as close to site average temperature as you can get)
bdot = .22 # average accumulation rate in m w.e. yr-1 -- make sure you convert your units if not originally in these units
rhos = 300 # surface density in kg m-3
z = np.arange(0.01,100,0.01) # depth array


# =============================================================================
# Below here is the model
# =============================================================================
Tavg = Tavg + 273.15    # Convert input temperature in celsius to Kelvin
bdot = bdot*rhoi    # Convert accumulation into units of m ice eq.

# Calculate the empirically fit slopes
c0=11*(bdot/rhow)*np.exp(-10160./(R*Tavg))
c1=575*np.sqrt(bdot/rhow)*np.exp(-21400./(R*Tavg))

# Calculate the two k constants
k0 = c0/bdot
k1 = c1/bdot

# Find the critical depth at which rho = rhoc
zc=(np.log(rhoc/(rhoi-rhoc)) - np.log(rhos/(rhoi-rhos)))/(k0*rhoi)

ix1 = z[z <= zc]                # find the z's above and below zc
ix2 = z[z > zc]
upix = np.searchsorted(z, ix1)  # indices above zc
dnix = np.searchsorted(z, ix2)  # indices below zc

q = np.zeros(len(z)) # pre-allocate some space for q
q[dnix] = np.exp(k1*rhoi*(z[dnix]-zc) + np.log(rhoc/(rhoi-rhoc)))
q[upix] = np.exp(k0*rhoi*z[upix] + np.log(rhos/(rhoi-rhos)))


# Here is where we produce the density profile!
rho = (q*rhoi)/(1 + q)


# create a porosity profile
#phi = 1-rho/rhoi

# =============================================================================
# Just make a basic plot to show the density profile
# =============================================================================
fig, ax = plt.subplots(figsize=(8,14))
ax.plot(rho, z, linewidth = 5, color = 'k')
ax.set_ylim(0, 100)
ax.invert_yaxis()
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.set_ylabel('Depth (m)', labelpad = 15, rotation = 270)
ax.set_xlabel('Density (kg m$^\mathrm{-3}$)')
# ax2 = ax.twiny()
# ax2.plot(phi, z, linewidth = 5, color = 'darkblue')
# ax2.set_xlabel('Porosity [-]', labelpad = 15)
# ax2.tick_params(axis='x', colors='darkblue') 
# ax2.xaxis.label.set_color('darkblue')
fig.tight_layout()
#ax2.tick_params(axis='x', which='major', pad=20)
#plt.savefig('density_porosity.png', transparent=True)