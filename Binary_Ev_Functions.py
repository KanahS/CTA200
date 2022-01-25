#!/usr/bin/env python
# coding: utf-8

# Functions to Plot Evolution of Binary: Eccentricity and Semi-Major Axis
# 
# only actually need to plot these values based on data from second binary as the primary is the center
# of the simulation
# 

# In[ ]:


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import rebound as rb
import reboundx as rx


# In[1]:


def binary_semi(archive, extras):
    sim_arc = rx.SimulationArchive(archive, rebxfilename=extras)
    
    x_arr = []
    y_arr = []
    
    for snap in range(len(sim_arc)): 
        base = sim_arc[snap][0].particles[1]
        orb_element = base.a
        time = (sim_arc[snap][0].t)
        
        y_arr.append(orb_element)
        x_arr.append(time)
        
    return x_arr, y_arr


# In[2]:


def binary_ecc(archive, extras):
    sim_arc = rx.SimulationArchive(archive, rebxfilename=extras)
    
    x_arr = []
    y_arr = []
    
    for snap in range(len(sim_arc)): 
        base = sim_arc[snap][0].particles[1]
        orb_element = base.e
        time = (sim_arc[snap][0].t)
        
        y_arr.append(orb_element)
        x_arr.append(time)

    return x_arr, y_arr


# In[ ]:




