from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.structure import Structure
from pyrasorbitalkit import OrbFs
from pyrasorbitalkit import wrp
from pyrasorbitalkit import ClcFs
import numpy as np

#- - - Computed the NOF-DFT energy for a given case
#- - - QCHEM KEYWORD VALUES FOR RASCI CALCULATION - - -
bas = 'cc-pVDZ'#'6-31g'
act = 8
elec = 6
occ = 4
spin_mult = 1
roots = 1
verbose=2
#- - - Molecule structure - - -
symbols = ['N', 'N']
chg = 0
mul = 1
coordinates = [[0.000, 0.000, 0.000],
               [0.000, 0.000, 1.100]]

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
#- - - Extract information for calculation of natural orbitals - - -
ovlp, mos, dm = wrp.rasci_ee_print(ee, False, False, False)

#- - - Compute Natural orbitals and their occupations - - -
nos, occs = OrbFs.Build_NOs(dm, mos, ovlp)
print('Natural Orbital occupations :', occs)

#- - - Print occupations in the active space - - -
occ_ras2, elec_ras2, occ_larger = OrbFs.RAS2_occupations(occs, occ, elec, act)
print('\nActive space NO occupations')
print('Occupations :', occs)
print('occ_ras2 = ', occ_ras2)
print('elec_ras2 = ', elec_ras2)
print('occ_low = ', occ_larger)


print('\nComputing RAS energy with Natural orbitals')
eNOF = ClcFs.RAS_NO(molecule, nos, bas, act, elec, occ,
                     spin_mult, roots, verbose)

print('\nRAS-NOF1 DFT energies : ', eNOF[0])

#- - - Compare with RASCI(h,p) total energy - - -
pars_rasci = parser_rasci(rascidata)
rasci_tot_e = pars_rasci['excited_states'][0]['total_energy']
print('Total RASCI(h,p) energy : ', rasci_tot_e)