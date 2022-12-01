#!/usr/bin/env python
# coding: utf-8

import reboundx as rx
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import rebound as rb


tup_num = 5 # change back to 50
e_b = np.linspace(0, 0.8, tup_num)
a_p = np.linspace(1, 5, tup_num)
Np = 15

Qex = [np.inf]
tidal_lag = [1,2,3,4]
for x in tidal_lag:
    Q = 10**x
    Qex.append(Q)

tup_list = []

for Q in Qex:
    for e in e_b:
        for a in a_p:
            tup_list.append((Q,e,a,Np))

Nq = len(Qex)
Ne = len(e_b)
Na = len(a_p)

def survival(initial):
    Q, eb, ap, Np = initial[0], initial[1], initial[2], initial[3]
    sim = rb.Simulation()
    sim.integrator = "whfast"

    mu = 0.5
    m1 = 1
    m2 = abs((m1*mu)/(1-mu))

    sim.add(m=m1, hash="Binary 1")
    sim.add(m=m2, a=1, e= eb, hash="Binary 2")

    R_star = 0.1*(1-0.5)**(1/3)
    sim.particles[0].r = R_star
    sim.particles[1].r = R_star

    #initializing Np massless planets
    for i in range(Np):
        f_plan = np.random.rand()*2*np.pi
        sim.add(m=0, a= ap, e=0, f= f_plan)

    #array to keep track of survival times
    sim.move_to_com()

    # Adding Tidal Elements
    rebx = rx.Extras(sim)
    tides = rebx.load_force("tides_constant_time_lag")
    rebx.add_force(tides)

    ps = sim.particles

    k2 = ps[0].params["tctl_k2"] = 0.035
    nb = ps[1].n


    if type(Q) != int:
        tau = ps[0].params["tctl_tau"] = 0
    else:
        tau = ps[0].params["tctl_tau"] = 3/(2*Q*k2*nb)
    
    # integrate
    N_times = int(10000)
    N_orbit = (1e4)*2*np.pi
    times = np.linspace(0,N_orbit,N_times)
    
    #array for survival times
    surv = np.zeros(Np)

    for i, time in enumerate(times):

        nb = ps[1].n
        r_eb = ps[1].e

        N_re = (1+(15./2.)*r_eb**2+(45./8.)*r_eb**4+(5./16.)*r_eb**6)/(1-r_eb**2)**6
        Omega_re = (1+3*r_eb**2+(3./8.)*r_eb**4)/(1-r_eb**2)**(9./2.)
        ps[0].params["Omega"] = N_re/Omega_re*nb

        sim.integrate(time, exact_finish_time=0)

        for num in reversed(range(2, sim.N)):
            p = sim.particles[num]
            o = p.calculate_orbit()
            p0 = sim.particles[0]
            p1 = sim.particles[1]
            d = np.sqrt(p.x**2 + p.y**2)
            d0 = np.sqrt((p.x - p0.x)**2 + (p.y - p0.y)**2)
            d1 = np.sqrt((p.x - p1.x)**2 + (p.y - p1.y)**2)

            if d > 20 or d0 < 0.25 or d1 < 0.25 or o.e > 1:
                surv[num-2] = time
                print(f'removing planet {num}') #
                sim.remove(num)

            if sim.N==2:
                break
    surv[(surv==0)] = time
    print(f'simulation finished, {len(sim.particles)-2} planets remaining')
    return np.mean(surv)

pool = rb.InterruptiblePool(processes=32)
mapping = pool.map(func= survival, iterable= tup_list)

def colour_plot(mapp):
    """
    Plots the colourmap from the small scale simulation
    """
    fig, ax = plt.subplots(1, Nq, figsize=(20,5), constrained_layout=True)
    ax = ax.ravel()
    
    SurvTimeAll = np.reshape(mapp, [Nq,Ne,Na])
    SurvTimeArr = [SurvTimeAll[i,:,:] for i in range(Nq)]
    
    for i in range(Nq):
        pcm = ax[i].pcolormesh(e_b, a_p, SurvTimeArr[i].T, shading='auto')
    
        a_b = 2.278 + 3.824*e_b - 1.71*(e_b**2)
        a_c = 1.6 + 5.1*e_b + (- 2.22*(e_b**2)) + 4.12*0.5 + (- 4.27*e_b*0.5) + (- 5.09*(0.5**2)) + 4.61*(e_b**2)*(0.5**2)
    
        ax[i].plot(e_b, a_c, color='lightsteelblue')
        ax[i].scatter(e_b, a_c, color='lightsteelblue')

        ax[i].plot(e_b, a_b, color='olive')
        ax[i].scatter(e_b, a_b, color='olive')
    
        ax[i].set_title('Q={:.1e}'.format(Qex[i]))
        ax[i].set_xlabel('Binary Eccentricity (e)')
        ax[i].set_ylabel('Planetary Semi-Major Axis (a)')
        #ax[i].set_xlim(0.0,0.8)
        #ax[i].set_ylim(1,5)

    plt.colorbar(pcm, location='right',label='Test Particle Survival Times')

colour_plot(mapping)




