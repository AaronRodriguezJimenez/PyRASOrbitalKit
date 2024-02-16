from pyqchem.structure import Structure
from pyqchem.parsers.parser_rasci import parser_rasci
from pyrasorbitalkit import ClcFs
from pyrasorbitalkit import OrbFs

"""
 In this example, we will perform a RASCI calculation with NOs followed by 
 a standard RASCI calculation and print the information.
"""

#- - - QCHEM KEYWORD VALUES FOR RASCI CALCULATION - - -
bas = 'cc-pVDZ'
act = 5
elec = 7
occ = 0
spin_mult = 4
roots = 3
verbose=2
#- - - Molecule structure - - -
symbols = ['N']
chg = 0
mul = 4
coordinates = [[0.000, 0.000, 0.000]]

molecule = Structure(coordinates=coordinates,
                     symbols=symbols,
                     charge=chg,
                     multiplicity=mul)

#- - -Perform RASCI QCHEM calculation and get information - - -
rascidata, ee = ClcFs.RASCICalculation(molecule, bas, act, elec, occ,
                                       spin_mult, roots, verbose)

print(rascidata)
print('\nInitial Active Space')
print('Active electrons : ', elec)
print('Active orbitals  : ', act)
print('Occ orbitals     : ', occ)

#- - - Read natural orbials and occupations from ee
nos = ee['nato_coefficients']
nos['qchem_order'] = ee['coefficients']['qchem_order']
occs = ee['nato_occupancies']['alpha']
# Get first NO occupation outside RAS2
occ_ras2, elec_ras2, occ_lower = OrbFs.RAS2_occupations(occs, occ, elec, act)
print('Occupations :', occs)
print('occ_ras2 = ', occ_ras2)
print('elec_ras2 = ', elec_ras2)
print('occ_low = ', occ_lower)

#- - - Perform RASCI calculation with NOs - - -
data, eerasnof, Elarger = ClcFs.RAS_NO(molecule, nos, bas, act, elec, occ,
                                       spin_mult, roots, verbose)
print('\nEnergy RAS with NOs : ', Elarger[0])

#- - - Compare with RASCI(h,p) total energy - - -
pars_rasci = parser_rasci(rascidata)
rasci_tot_e = pars_rasci['excited_states'][0]['total_energy']
print('\nTotal RASCI(h,p) energy : ', rasci_tot_e)
print('\nTotal E_RAS_NOs_large :', Elarger[0][0])