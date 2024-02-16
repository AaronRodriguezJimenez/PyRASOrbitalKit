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
print('occupations :', occs)
#- - - Redefine RAS2 based on occupations? - - -
occ_lev = 0.01
act_new, elec_new, occ_new, occ_lower = OrbFs.NewRAS2(occs, occ_lev)
print('\nRedefined Active Space')
print('New active electrons : ', elec_new)
print('New active orbitals  : ', act_new)
print('New occ orbitals     : ', occ_new)
print('Lowest NO occ in RAS2: ', occ_lower)

#- - - Avoid redefinition of RAS2 but cut orbital space there
# Get first NO occupation outside RAS2
occ_ras2, elec_ras2, occ_larger = OrbFs.RAS2_occupations(occs, occ, elec, act)
print('Occupations :', occs)
print('occ_ras2 = ', occ_ras2)
print('elec_ras2 = ', elec_ras2)
print('occ_low = ', occ_larger)

#- - - Compute DFT energy for a given threshold - - -
#eNOF is a list with DFT values for each root
print('\nComputing NOF1-DFT energy for lower occs than RAS2')
eNOF = ClcFs.NOFDFT2_savin(molecule, occ_lower, bas, act, elec, occ,
                     spin_mult, roots, verbose)

print('\nRAS-NOF1 DFT energies : ', eNOF[0])

#- - - Compare with RASCI(h,p) total energy - - -
pars_rasci = parser_rasci(rascidata)
rasci_tot_e = pars_rasci['excited_states'][0]['total_energy']
print('Total RASCI(h,p) energy : ', rasci_tot_e)