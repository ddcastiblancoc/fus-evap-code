#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# >> Python libraries used in this code. fisbar needs to be installed manually; see: gitlab.in2p3.fr/gregoire.henning/fisbar-python/-/tree/v001.
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.misc import derivative
from fisbar.fisbar import fisbar
import random as rn
from joblib import Parallel, delayed
import multiprocessing
from collections import Counter
import time

# >> Data available of all isotopes (mass number, atomic number, element name, mass excess,...) in NuDat database, NNDC; see: www.nndc.bnl.gov/nudat3/indx_sigma.jsp.
isotopes_chart = pd.read_csv('isotopes_data_NuDat.txt', delimiter = "\t")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# >> atomic_number returns the atomic number given the element name as a string.
def atomic_number(element):
    i = isotopes_chart[isotopes_chart.iloc[:,1].str.strip()==element]
    if i.shape[0]!=0:
        return int(i.iloc[0,2])
        
# >> element_number returns the element name as a string given the atomic number Z.
def element_symbol(Z):
    i = isotopes_chart[isotopes_chart.iloc[:,2]==Z]
    if i.shape[0]!=0:
        return i.iloc[0,1].strip()
        
# >> ME returns the mass excess of nucleus A,Z.
def ME(A,Z):
    i = isotopes_chart[(isotopes_chart.iloc[:,0]==A) & (isotopes_chart.iloc[:,2]==Z)]
    if i.shape[0]!=0:
        return float(i.iloc[0,6])
        
# >> nuclear_mass returns mass of nucleus A,Z from mass excess value.
def nuclear_mass(A,Z):
    return ME(A,Z)+A*931.5
    
# >> kinematic_quantities returns some physical quantities (CM energy, Q value, excitation energy,...) related to the formation of compound nucleus given projectile 
# >> and target nucleus.
def kinematic_quantities(A_proj, A_targ, Z_proj, Z_targ, E_inc_lab):
    m_proj = nuclear_mass(A_proj, Z_proj)
    m_targ = nuclear_mass(A_targ, Z_targ)
    m_comp = nuclear_mass(A_proj+A_targ, Z_proj+Z_targ)
    
    E_inc_CM = (m_targ*E_inc_lab)/(m_targ+m_proj)
    r_proj = 1.2*np.cbrt(A_proj)
    r_targ = 1.2*np.cbrt(A_targ)
    coulomb_barrier = 1.44*Z_targ*Z_proj/(r_proj+r_targ)
    Q_compound = m_proj+m_targ-m_comp
    excit_energy = Q_compound+E_inc_CM
    v_proj_lab_c = np.sqrt(2*E_inc_lab/m_proj)
    v_comp_lab_c = m_proj*v_proj_lab_c/m_comp
    
    return E_inc_CM, coulomb_barrier, Q_compound, excit_energy, v_proj_lab_c, v_comp_lab_c
    
# >> fusion_potential returns Bass + Coulomb potential between projectile and target at a given distance r.
def fusion_potential(r, A_proj, A_targ, Z_proj, Z_targ):
    r1 = 1.16*np.cbrt(A_proj)
    r2 = 1.16*np.cbrt(A_targ)
    R1 = r1-1.6124/r1
    R2 = r2-1.6124/r2
    s = r-R1-R2
    g = 1/(0.03*np.exp(s/3.30)+0.0061*np.exp(s/0.65))
    
    bass_potential = -(R1*R2/(R1+R2))*g
    coulomb_potential = 1.44*Z_targ*Z_proj/r
    
    return bass_potential + coulomb_potential
    
# >> fusion_barrier returns position and height of Bass barrier given projectile and target nucleus.
def fusion_barrier(A_proj, A_targ, Z_proj, Z_targ):    
    r_proj = 1.2*np.cbrt(A_proj)
    r_targ = 1.2*np.cbrt(A_targ)
    r_barr = minimize(lambda r: -fusion_potential(r, A_proj, A_targ, Z_proj, Z_targ), r_proj+r_targ)
    
    return r_barr.x[0], fusion_potential(r_barr.x[0], A_proj, A_targ, Z_proj, Z_targ)

# >> classical_cs returns classical fusion cross section given projectile, target nucleus and CM incident energy; includes correction by l_max from fisbar.
def classical_cs(E, A_proj, Z_proj, A_targ, Z_targ):
    R, V = fusion_barrier(A_proj, A_targ, Z_proj, Z_targ)
    m_proj = ME(A_proj,Z_proj)+A_proj*931.5
    m_targ = ME(A_targ,Z_targ)+A_targ*931.5
    reduced_m = m_targ*m_proj/(m_targ+m_proj)
    
    fisdata = fisbar(Z_proj+Z_targ, A_proj+A_targ, 0)
    l_crit = fisdata['elmax']
    x = 197.3**2*l_crit**2/(2*reduced_m)
    E_crit = V + x/R**2
    
    if E<V:
        return 0
    else:
        if E>E_crit and fisdata['success']==True:
            return 10*np.pi*x/E
        else:
            return 10*np.pi*R**2*(1-(V/E))

# >> Wong_cs returns fusion cross section by Wong formula given projectile, target nucleus and CM incident energy; includes correction by l_max from fisbar.
def Wong_cs(E, A_proj, Z_proj, A_targ, Z_targ):
    R, V = fusion_barrier(A_proj, A_targ, Z_proj, Z_targ)
    m_proj = ME(A_proj,Z_proj)+A_proj*931.5
    m_targ = ME(A_targ,Z_targ)+A_targ*931.5
    reduced_m = m_targ*m_proj/(m_targ+m_proj)
    
    V_tot = lambda r: fusion_potential(r, A_proj, A_targ, Z_proj, Z_targ)
    dV_tot = lambda r: derivative(V_tot, r, dx=1e-3)
    ddV_tot = lambda r: derivative(dV_tot, r, dx=1e-3)
    hbar_omega = 197.3*np.sqrt(-ddV_tot(R)/reduced_m)
    
    fisdata = fisbar(Z_proj+Z_targ, A_proj+A_targ, 0)
    l_crit = fisdata['elmax']
    x = 197.3**2*l_crit**2/(2*reduced_m)
    E_crit = V + x/R**2
    
    if E>E_crit and fisdata['success']==True:
        return 10*np.pi*x/E
    else:
        return (10*hbar_omega*R**2/(2*E))*np.log(1+np.exp((2*np.pi/hbar_omega)*(E-V)))

# >> level_density returns level density (Gilbert-Cameron) for a nucleus A,Z and excitation energy E.
def level_density(A,Z,E):
    a = A/8 
    return (np.sqrt(np.pi)*np.exp(2*np.sqrt(a*E)))/(12*a**0.25*E**1.25)

# >> Weisskopf returns probability distribution for ejectile energy (neutron, proton and alpha) and the corresponding energy range, given a compound nucleus A,Z
# >> and a excitation energy E.
def Weisskopf(A, Z, E, ejectile):    
    if ejectile=='neutron':
        A_ejec = 1
        Z_ejec = 0
        s_j = 0.5
    elif ejectile=='proton':
        A_ejec = 1
        Z_ejec = 1
        s_j = 0.5
    elif ejectile=='alpha':
        A_ejec = 4
        Z_ejec = 2
        s_j = 0
    A_resid = A-A_ejec
    Z_resid = Z-Z_ejec
    
    ME_ejec = ME(A_ejec,Z_ejec)
    Q = ME(A,Z)-ME_ejec-ME(A_resid,Z_resid)
    g_j = (2*s_j+1)*(ME_ejec+A_ejec*931.5)/(np.pi**2*197.3**2)
    R_g = 1.2*(np.cbrt(A_ejec)+np.cbrt(A_resid))
    sigma_g = np.pi*R_g**2
    V_coulomb = 1.44*Z_ejec*Z_resid/R_g
    
    if ejectile=='neutron':
        beta = 0
    else:
        beta = -V_coulomb
    
    if V_coulomb<E+Q:
        e = np.linspace(V_coulomb,E+Q,100)
        e = e[:-1]
    else:
        e = None
    
    P = g_j*sigma_g*(e+beta)*level_density(A_resid,Z_resid,E+Q-e)/level_density(A,Z,E)
    return e, P

# >> BohrWheeler returns probability distribution for deformation energy (fission), given a compound nucleus A,Z and a excitation energy E.
def BohrWheeler(A, Z, E):
    fisdata = fisbar(Z, A, 0)
    if fisdata['success']==True:
        B_fis = fisdata['bfis']
    else:
        B_fis = None
    
    if 0<E-B_fis:
        e = np.linspace(0,E-B_fis,100)
        e = e[:-1]
    else:
        e = None
    
    P = level_density(A,Z,E-B_fis-e)/(2*np.pi*level_density(A,Z,E)*(1+np.exp(2*np.pi*(B_fis-e))))
    return e, P

# >> branching_ratio returns decay width for both evaporation (Weisskopf) and fission (BohrWheeler) integrating the probability distribution; also returns
# >> normalized probability distribution for energy.
def branching_ratio(A, Z, E, ejectile):
    if ejectile=='fission':
        e, P = BohrWheeler(A, Z, E)
    else:
        e, P = Weisskopf(A, Z, E, ejectile)
    dx = e[1]-e[0]
    G = sum(P*dx)
    return G, e, P/G

# >> energy_ejectile returns a random value of ejectile energy using Monte Carlo with the energy probability distribution (normalized).
def energy_ejectile(x, pdf):
    dx = x[1]-x[0]
    F = np.cumsum(pdf*dx)
    
    r = rn.uniform(F[0],F[-1])
    for j in range(1,F.size):
        if r >= F[j-1] and r < F[j]:
            return x[j]
            break    

# >> decay returns a residue nucleus after choosing a decay randomly between neutron, proton, alpha emission and fission with Monte Carlo. If fission is chosen, 
# >> it returns just a 'fission' string. If a decay width can't be calculated, it's forced to be 0. If all of them are 0, none of the decay modes can occur: the residue
# >> nuclei doesn't change
def decay(A, Z, E):
    try:
        G_n, e_n, P_n = branching_ratio(A, Z, E, 'neutron')
    except:
        G_n = 0
    try:
        G_p, e_p, P_p = branching_ratio(A, Z, E, 'proton')
    except:
        G_p = 0
    try:
        G_a, e_a, P_a = branching_ratio(A, Z, E, 'alpha')
    except:
        G_a = 0
    try:
        G_f, e_f, P_f = branching_ratio(A, Z, E, 'fission')
    except:
        G_f = 0
    
    G = G_n+G_p+G_a+G_f
    if G==0:
        return A, Z, E
    
    r = rn.random()
    if r<G_n/G:
        Q = ME(A,Z)-ME(1,0)-ME(A-1,Z)
        A=A-1
        Z=Z
        E = E+Q-energy_ejectile(e_n, P_n)
        return A, Z, E
    
    elif G_n/G<r<(G_n+G_p)/G:
        Q = ME(A,Z)-ME(1,1)-ME(A-1,Z-1)
        A=A-1
        Z=Z-1
        E = E+Q-energy_ejectile(e_p, P_p)
        return A, Z, E
    
    elif (G_n+G_p)/G<r<(G_n+G_p+G_a)/G:
        Q = ME(A,Z)-ME(4,2)-ME(A-4,Z-2)
        A=A-4
        Z=Z-2
        E = E+Q-energy_ejectile(e_a, P_a)
        return A, Z, E
    
    elif (G_n+G_p+G_a)/G<r<1:
        return 'fission'

# >> cascade returns the final nucleus from a cascade of decays. If fission is selected, the cycle stops. If none of the modes is available for decay (when the residue
# >> doesn't change after a decay), the cycle stops.
def cascade(A, Z, E):
    while True:
        d = decay(A, Z, E)
        if d=='fission':
            return 'fission'
            break
        else:
            A1, Z1, E1 = d
            if E1 == E:
                return (A, Z)
                break
            A, Z, E = A1, Z1, E1
            
# >> cascade_counts returns the list of all the exit channels with each corresponding number of events: with Counter, all the events of all cascades are organized
# >> counting how many common events are found of each nucleus (or fission event).
def cascade_counts(A, Z, E, N_casc):
    inputs = range(N_casc)
    processInput = lambda i: cascade(A, Z, E)

    num_cores = multiprocessing.cpu_count()
    results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)
    count = Counter(results)
    
    return count

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------