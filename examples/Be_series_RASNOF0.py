from pyrasorbitalkit.parser_rasci_rasnof0 import parser_rasci
from pyqchem.structure import Structure
from pyrasorbitalkit import ClcFs
import numpy as np

"""
 This script performs calcualtions within the RASNOF0 approximation,
 check availability of the corresponding qchem keywords before running the script.
"""

#* * * CALCULATIONS FOR HE-SERIES * * *
E_RAS = []
E_RAS_NOs = []
fci_energy = []

#- - - QCHEM KEYWORD VALUES FOR RASCI CALCULATION - - -
bas = 'cc-pvqz'#'sto-3g'
act = 5
elec = 4
occ = 0
spin_mult = 0
roots = 1
verbose=2
occthrs=0.001
version=3
functional='VWN1_NOF3'

#- - - Series symbols - - -
be_series = ['Li','Be','B','C','N','O','F','Ne']
chg = -1
mul = 1
for atom in be_series:
    #Build QC input
    molecule = Structure(coordinates=[[0.000,0.000,0.000]],
                         symbols=[atom],
                         charge=chg,
                         multiplicity=mul)

    #- - -Perform RASCI QCHEM calculation and get information - - -
    print('\n# # # # # CALCULATION FOR %s+%d ATOM # # # # #' %(atom, chg))
    #Threshold are checked previously, work for cc-pVQZ basis set
    if atom == 'Li':
        occthrs = 0.025
    if atom == 'Be':
        occthrs = 0.0097
    if atom == 'B':
        occthrs = 0.0016
    if atom == 'C':
        occthrs = 0.0004
    if atom == 'N':
        occthrs = 0.00017
    if atom == 'O':
        occthrs = 0.00011
    if atom == 'F':
        occthrs = 0.000095
    if atom == 'Ne':
        occthrs = 0.000095
    #- Calculation for RASNOF1,2 versions
    #rascidata, ERAS, ERASNOF = ClcFs.RASNOF0(molecule, occthrs, version, functional, bas,
    #                                              act, elec, occ, spin_mult, roots, verbose)
    #- Calculation for RASNOF3 version
    rascidata, ERAS, ERASNOF = ClcFs.RASNOF30_sum(molecule, occthrs, version, functional, bas,
                                                  act, elec, occ, spin_mult, roots, verbose)
    print(rascidata)
    print(ERAS[0][0])
    E_RAS.append(ERAS[0][0])
    print(ERASNOF[0][0])
    E_RAS_NOs.append(ERASNOF[0][0])

    #- - -Perform a FCI/cc-pVDZ calculation for atom
    if atom=='He':
        fcidata, eefci = ClcFs.RASCICalculation(molecule, bas, 5, elec,occ, spin_mult, roots, verbose)
        fcidata = parser_rasci(fcidata)
        print(fcidata)
        fci_energy.append([state['total_energy'] for state in fcidata['excited_states']])

    else:
        fcidata, eefci = ClcFs.RASCICalculation(molecule, bas, 14, elec,occ, spin_mult, roots, verbose)
        fcidata = parser_rasci(fcidata)
        fci_energy.append([state['total_energy'] for state in fcidata['excited_states']])


    print('\n# # # # # END OF CALCULATION for %s+%d ATOM # # # # #' %(atom, chg))
    print('SUM THRESHOLD :', occthrs)
    chg+=1

# transform energy list in array
E_RAS = np.array(E_RAS)
E_RAS_NOs = np.array(E_RAS_NOs)
print(fci_energy)

# print energies
print('---------- R E S U L T S -----------')
if version > 2:
    print('Sum Threshold = ', occthrs)
    print('\n# Atom  RASCI(h,p)   EFCI      RAS-NOF0')
    print('------------------------------------------')
    for d, e, efci, enof in zip(be_series, E_RAS, fci_energy, E_RAS_NOs):
        print('{:s} : {:10.5f} : {:10.5f} : {:10.5f}'.format(d, e, efci[0], enof))

    for enof in E_RAS_NOs:
        print('{:10.5f}'.format(enof))
else:
    print('Occupation Threshold = ', occthrs)
    print('\n# Atom  RASCI(h,p)   EFCI      RAS-NOF0')
    print('------------------------------------------')
    for d, e, efci, enof in zip(be_series, E_RAS, fci_energy, E_RAS_NOs):
        print('{:s} : {:10.5f} : {:10.5f} : {:10.5f}'.format(d, e, efci[0], enof))

    for enof in E_RAS_NOs:
        print('{:10.5f}'.format(enof))