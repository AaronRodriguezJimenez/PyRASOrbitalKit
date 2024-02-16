import numpy as np
from pyrasorbitalkit import ClcFs
from pyqchem.structure import Structure
from pyrasorbitalkit import OrbFs as oei

"""
    The nu values are the threshold for the natural orbitals occupation.
    The nu values are used to compute the RASCI energy with the NOs as initial guess.
    Natural orbital occupation numbers are computed from the 1-PDM, the MOs and S.
"""

#- - - QCHEM KEYWORD VALUES FOR RASCI CALCULATION - - -
bas = 'cc-pvdz'#'sto-3g'
act = 6
elec = 6
occ = 4
spin_mult = 1
roots = 1
verbose=2
E_RAS = []
E_RAS_NOs = []
E_DFT = []

# - - - Structures - - -
ethylene_0deg_xyz = [[ 1.3131719011, 0.0000000000,  0.9253935281],
           [ 1.3131667074, 0.0000000000, -0.9253979612],
           [-1.3131719011, 0.0000000000,  0.9253935281],
           [-1.3131667074, 0.0000000000, -0.9253979612],
           [ 0.7753000000, 0.0000000000, -0.0000009008],
           [-0.7753000000, 0.0000000000, -0.0000009008]]
ethylene_symbols = ['H', 'H', 'H', 'H', 'C', 'C']

be_xyz = [[0.0000000000, 0.0000000000, 0.0000000000]]
be_symbols = ['Be']

H2_xyz = [[0.0000000000, 0.0000000000, 0.3700000000],
          [0.0000000000, 0.0000000000,-0.3700000000]]
H2_symbols = ['H', 'H']

N2_xyz = [[0.0000000000, 0.0000000000, 0.5500000000],
            [0.0000000000, 0.0000000000,-0.5500000000]]
N2_symbols = ['N', 'N']

#- - - - - -
molecule = Structure(coordinates=N2_xyz,
                     symbols=N2_symbols,
                     charge=0,
                     multiplicity=3)

rascidata, ee = ClcFs.RASCICalculation(molecule, bas, act, elec, occ, spin_mult, roots, verbose)
S = np.array(ee['overlap'])
DM1 = np.array(ee['total_scf_density_multi'][0]) #for state 1 for more states needs changes
C = np.array(ee['coefficients']['alpha'])
nos, noons = oei.Build_NOs(DM1, C, S)
#Electronic structure NOONs
#nos = ee['nato_coefficients']
#nos['qchem_order'] = ee['coefficients']['qchem_order']
#noons = ee['nato_occupancies']['alpha']
#print(occs)

print('- - - - - - - - - - - - - - ')
state = 0 #state label
nu_list = []
total_orbs = len(noons)
print('Total orbitals : ', total_orbs)
occ, elec, occ_low = oei.RAS2_occupations(noons, occ, elec, act)
print('First occupation at RAS3 : ', occ_low)

for i, nu in enumerate(reversed(noons)):
    nu_list.append(nu)
    print('Frozen Virtuals : ', i+1, nu)
    if (nu > occ_low):
        print('nu Max = ', nu)
        break


#- - - print & plot results
E_RAS = np.array(E_RAS)
E_RAS_NOs = np.array(E_RAS_NOs)
print('\n#occupation  RASCI(6,6)   RAS-NO>nu')
print('----------------------------')
for d, e, elarge in zip(nu_list, E_RAS, E_RAS_NOs):
    print('{:4.1f} : {:10.5f} {:10.5f}'.format(d, e, elarge))