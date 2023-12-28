#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# The purpose of this Python code is calculate cross sections of each exit channel of a fusion-evaporation reaction, given a projectile and a target.
# This code uses Bass model and Wong formula (or classical) to calculate fusion cross section, Monte Carlo to select a decay mode randomly, Weisskopf-Ewing formula 
# to calculate evaporation decay widths and Bohr-Wheeler model to calculate fission decay width.
# It receives as input: mass number of projectile, atomic number of projectile, mass number of target, atomic number of target, laboratory energy of projectile, and
# number of cascades.
# Compiling this code returns (in console): (1) data about projectile, target and compound nucleus, (2) physical quantities related to the formation of compound 
# nucleus, (3) list of found residues from compound evaporation (and each evaporation cross section).

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# >> Importation of routines and libraries used in this code
from routines_fus_evap_code import *

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# >> Input information
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

print('Insert A_proj, sym_proj, A_targ, sym_targ, E_inc_lab, N_casc, separated by spaces.')
print('For example, for a 58-Ni beam on 50-Cr at 250 MeV and 1000 cascades, insert: 58 Ni 50 Cr 250 1000.')
print('I recommend using 1000 cascades (it\'s going to take 1-2 minutes), I don\'t recommend more than 10000 cascades, too much time!')

A_proj, sym_proj, A_targ, sym_targ, E_inc_lab, N_casc = input().split(' ')
A_proj, A_targ, E_inc_lab, N_casc = int(A_proj), int(A_targ), float(E_inc_lab), int(N_casc)
Z_proj, Z_targ = atomic_number(sym_proj), atomic_number(sym_targ)

print('\nCalculating...\n')

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# >> Calculation of physical quantities and list of exit channels
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

t_i=time.time()

E_inc_CM, coulomb_barrier, Q_compound, excit_energy, v_proj_lab_c, v_comp_lab_c = kinematic_quantities(A_proj, A_targ, Z_proj, Z_targ, E_inc_lab)
R_barrier, V_barrier = fusion_barrier(A_proj, A_targ, Z_proj, Z_targ)

fusion_cs = Wong_cs(E_inc_CM, A_proj, Z_proj, A_targ, Z_targ)

count = cascade_counts(A_proj+A_targ, Z_proj+Z_targ, excit_energy, N_casc)

exit_channels = [key for key in count]
A_res = np.array([item[0] for item in exit_channels if item!='fission'])
Z_res = np.array([item[1] for item in exit_channels if item!='fission'])
N_res = A_res-Z_res
sym_res = np.vectorize(element_symbol)(Z_res)
events_res = np.array([count[x] for x in exit_channels if x!='fission'])
percent_res = 100*events_res/N_casc
cs_res = fusion_cs*percent_res/100

d = {'A':A_res, 'sym':sym_res, 'Z':Z_res, 'N':N_res, 'events':events_res, '%':percent_res, 'cross section (mb)':cs_res}
residual_nuclei = pd.DataFrame(data = d)
residual_nuclei = residual_nuclei.sort_values(by=['A','Z'], ascending=[True,True])
residual_nuclei.index = range(len(residual_nuclei))

if 'fission' in exit_channels:
    events_fis = count['fission']
    percent_fis = 100*events_fis/N_casc
    cs_fis = fusion_cs*percent_fis/100
    new_row = {'A':'fission', 'sym':'-', 'Z':'-', 'N':'-', 'events':events_fis, '%':percent_fis, 'cross section (mb)':cs_fis}
    residual_nuclei = pd.concat([residual_nuclei, pd.DataFrame([new_row])], ignore_index=True)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------
# >> Printing of results in console
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------

new_row = {'A':'total', 'sym':'-', 'Z':'-', 'N':'-', 'events':N_casc, '%':100, 'cross section (mb)':fusion_cs}
residual_nuclei = pd.concat([residual_nuclei, pd.DataFrame([new_row])], ignore_index=True)

output_1 = pd.DataFrame([
            [A_proj, sym_proj, Z_proj, A_proj-Z_proj, ME(A_proj,Z_proj)],
            [A_targ, sym_targ, Z_targ, A_targ-Z_targ, ME(A_targ,Z_targ)],
            [A_targ+A_proj, element_symbol(Z_targ+Z_proj), Z_targ+Z_proj, A_targ+A_proj-Z_targ-Z_proj, ME(A_targ+A_proj,Z_targ+Z_proj)]],
            columns = ['A', 'sym', 'Z', 'N', 'ME (MeV)'],
            index = ['Projectile', 'Target', 'Compound'])

output_2 = pd.DataFrame([
            [E_inc_lab], [E_inc_CM], [coulomb_barrier], [Q_compound], [excit_energy], [v_proj_lab_c*100], [v_comp_lab_c*100], [V_barrier], [R_barrier], [fusion_cs]],
            index = [
                    'Projectile energy in lab frame (MeV)',
                    'Center of mass projectile energy (MeV)',
                    'Projectile-target Coulomb barrier (MeV)',
                    'Q-value of compound reaction (MeV)',
                    'Excitation energy of compound nucleus (MeV)',
                    'Laboratory projectile nucleus velocity/c (%)',
                    'Laboratory compound nucleus velocity/c (%)',
                    'Bass barrier (nuclear+Coulomb) (MeV)',
                    'Bass barrier position (fm)',
                    'Compound formation cross section (mb)'
                    ])

print('>> Nuclei involved in fusion reaction')
print(output_1.to_string())
print('\n>> Compound formation related quantities')
print(output_2.to_string(header=False))
print('\n>> Found residues from compound evaporation')
print(residual_nuclei.to_string(index=False))

print('\nDone! It took ' + str(round(time.time()-t_i,3)) + ' seconds to generate the results.')

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------