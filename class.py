#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import rebound as rb

tup_num = 25
e_b = np.linspace(0, 0.7, tup_num)
a_p = np.linspace(1, 5, tup_num)
mass = np.linspace(0.1, 1.0, 10)
Np = 15
tup_list = []

for mu in mass:
    for e in e_b:
        for a in a_p:
            tup_list.append((mu,e,a,Np))
Nm = len(mass)
Ne = len(e_b)
Na = len(a_p)

def survival(initial):    
    mu, eb, ap, Np = initial[0], initial[1], initial[2], initial[3]

    sim = rb.Simulation()
    sim.integrator = "whfast"
    
    m1 = 1
    if mu==1:
        m2=m1
    else:
        m2 = abs((m1*mu)/(1-mu))
    
    sim.add(m=m1, hash="Binary 1")
    sim.add(m=m2, a=1, e= eb, hash="Binary 2")
    
    #initializing Np massless planets
    for p in range(Np):
        f_plan = np.random.rand()*2.*np.pi
        sim.add(m=0, a= ap, e=0, f= f_plan)
    
    #array to keep track of survival times
    sim.move_to_com()

    #directory_orbit = '/mnt/raid-cita/ksmith/ClassOrbParamsKC/'
    #filename_orbit = r"KCeb{:.3f}_ap{:.3f}_Np{:.1f}_tup{:.1f}_mu{:.1f}.bin".format(eb,ap,Np,tup_num,mu)
    #sim.automateSimulationArchive(directory_orbit+filename_orbit, interval=1e1, deletefile=True)

    
    #integrate
    N_times = int(100)
    N_orbit = (1e4)*2*np.pi
    times = np.linspace(0,N_orbit,N_times)
   # array for survival times
    surv = np.zeros(Np)

    for i, time in enumerate(times):
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
time_surv = np.reshape(mapping, [Nm,Ne,Na])

directory_surv = '/mnt/raid-cita/ksmith/ClassSurvTimes/'
txt_surv = f'map_tup{tup_num}plan{Np}_mu{mass[0]}_{mass[-1]}.txt'
npy_surv = f'map_tup{tup_num}plan{Np}_mu{mass[0]}_{mass[-1]}.npy'
bin_surv = f'map_tup{tup_num}plan{Np}_mu{mass[0]}_{mass[-1]}.bin'
np.savetxt(directory_surv+txt_surv, mapping)
np.savetxt(directory_surv+npy_surv, mapping)
np.savetxt(directory_surv+bin_surv, mapping)

figure_all = time_surv = np.reshape(mapping, [Nm,Ne,Na])
i = 0
figure = figure_all[i,:,:]

plt.pcolormesh(e_b, a_p, figure.T, shading='auto')
plt.title(f'Mean Survival Times (mu={mass[i]})')
plt.xlabel('Binary Eccentricity (e)')
plt.ylabel('Planetary Semi-Major Axis (a)')
plt.xlim(0.0,0.7)
plt.ylim(1,5)

a_c = 1.6 + 5.1*e_b - 2.22*(e_b**2) + 4.12*mass[i] - 4.27*e_b*mass[i] - 5.09*(mass[i]**2) + 4.61*(e_b**2)*(mass[i]**2)

a_b = 2.278 + 3.824*e_b - 1.71*(e_b**2)

plt.plot(e_b, a_c, color='lightsteelblue')
plt.scatter(e_b, a_c, color='lightsteelblue')
plt.plot(e_b, a_b, color='white')
plt.scatter(e_b, a_b, color='white')

plt.colorbar(label='Test Particle Survival Times')
plt.show()

directory_test = '/mnt/raid-cita/ksmith/'
completed = 'The simulation finished!'
name = 'DONE'
np.save(directory_test+name, completed)


