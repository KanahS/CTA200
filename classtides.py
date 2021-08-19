#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import rebound as rb
import reboundx as rx

tup_num = 7
e_b = np.linspace(0, 0.7, tup_num)
a_p = np.linspace(1, 5, tup_num)
Qex = []
for x in range(4,7):
    Q = 10**x
    Qex.append(Q)

Np = 5
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
    
    #initializing Np massless planets
    for p in range(Np):
        f_plan = np.random.rand()*2.*np.pi
        sim.add(m=0, a= ap, e=0, f= f_plan)
    
    #array to keep track of survival times
    sim.move_to_com()

    # Adding Tidal Elements 
    rebx = rx.Extras(sim)
    tides = rebx.load_force("tides_constant_time_lag")
    rebx.add_force(tides)
    
    ps = sim.particles
    
    k2 = ps[0].params["tctl_k1"] = 0.035
    nb = ps[1].n
        
    ps[0].params["tctl_tau"] = 3/(2*Q*k2*nb)
    tau = ps[0].params["tctl_tau"]

    directory_orbit = '/mnt/raid-cita/ksmith/ClassOrbParamsTidesKC/'
    filename_orbit = r"KCeb{:.3f}_ap{:.3f}_Np{:.1f}_tup{:.1f}_tau{:.4f}.bin".format(eb,ap,Np,tup_num,tau)
    sim.automateSimulationArchive(directory_orbit+filename_orbit, interval=1e1, deletefile=True)

    
    #integrate
    N_times = int(100)
    N_orbit = (1e4)*2*np.pi
    times = np.linspace(0,N_orbit,N_times)
   
    #array for survival times
    surv = np.zeros(Np)

    for i, time in enumerate(times):
   
        nb = ps[1].n
        r_eb = ps[1].e
        
        N_re = (1.+(15./2.)*r_eb**2+(45./8.)*r_eb**4+(5./16.)*r_eb**6)/(1-r_eb**2)**6
        Omega_re = (1+3*r_eb**2+(3./8.)*r_eb**4)/(1-r_eb**2)**(9./2.)
        ps[0].params["Omega"] = N_re/Omega_re*nb

        sim.integrate(time, exact_finish_time=0)

        for num in reversed(range(2, sim.N)):
            p = sim.particles[num]
            if (p.x**2 + p.y**2) > (100)**2:
                surv[num-2] = time
                print(f'removing planet {num}')
                sim.remove(num)

        if sim.N==2:
            break
    surv[(surv==0)] = time
    
    print(f'simulation finished, {len(sim.particles)-2} planets remaining')
    return np.mean(surv)
   
pool = rb.InterruptiblePool(processes=16)
mapping = pool.map(func= survival, iterable= tup_list)
time_surv = np.reshape(mapping, [Nq,Ne,Na])

directory_surv = '/mnt/raid-cita/ksmith/ClassSurvTimesTides/'
txt_surv = f'map_tup{tup_num}plan{Np}.txt'
npy_surv = f'map_tup{tup_num}plan{Np}.npy'
bin_surv = f'map_tup{tup_num}plan{Np}.bin'
np.savetxt(directory_surv+txt_surv, mapping)
np.savetxt(directory_surv+npy_surv, mapping)
np.savetxt(directory_surv+bin_surv, mapping)

fig, ax = plt.subplots(1, Nq, figsize=(20,5), constrained_layout=True)
ax = ax.ravel()
SurvTimeArr = [time_surv[i,:,:] for i in range(Nq)]

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
    ax[i].set_xlim(0.0,0.7)
    ax[i].set_ylim(1,5)

plt.colorbar(pcm, location='right',label='Test Particle Survival Times') 

# older plots

#figure_all = time_surv
#i = 0
#figure = figure_all[i,:,:]

#plt.pcolormesh(e_b, a_p, figure.T, shading='auto')
#plt.title(f'Mean Survival Times (Q={Qex[i]})')
#plt.xlabel('Binary Eccentricity (e)')
#plt.ylabel('Planetary Semi-Major Axis (a)')
#plt.xlim(0.0,0.7)
#plt.ylim(1,5)

#a_c = 1.6 + 5.1*e_b - 2.22*(e_b**2) + 4.12*0.5 - 4.27*e_b*0.5 - 5.09*(0.5**2) + 4.61*(e_b**2)*(0.5**2)
#a_b = 2.278 + 3.824*e_b - 1.71*(e_b**2)

#plt.plot(e_b, a_c, color='lightsteelblue')
#plt.scatter(e_b, a_c, color='lightsteelblue')
#plt.plot(e_b, a_b, color='white')
#plt.scatter(e_b, a_b, color='white')

#plt.colorbar(label='Test Particle Survival Times')
plt.show()

directory_test = '/mnt/raid-cita/ksmith/'
completed = 'The simulation finished!'
name = 'DONE'
np.save(directory_test+name, completed)
