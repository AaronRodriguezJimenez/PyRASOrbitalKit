from pyqchem.parsers.parser_rasci import parser_rasci
from pyqchem.structure import Structure
from pyrasorbitalkit import OrbFs as oei
from pyrasorbitalkit import wrp
from pyrasorbitalkit import ClcFs
import numpy as np

#- - - QCHEM KEYWORD VALUES FOR RASCI CALCULATION - - - 
bas = 'DZ*'#'sto-3g'
act = 6
elec = 6
occ = 4
spin_mult = 1
roots = 1
verbose=2
#- - - N-N separations - - -
short_range = np.arange(0.50, 0.90, 0.20)
bond_range= np.arange(1.0, 2.1, 0.1)
long_range = np.array([2.50, 3.00, 3.50, 4.00])
distances = np.concatenate((short_range, bond_range,long_range))

E_RAS = []
E_RAS_NOs = []
for b in distances:
    #- - - Molecule structure - - -
    symbols = ['N', 'N']
    chg = 0
    mul = 1
    coordinates = [[0.000, 0.000, 0.000],
                   [0.000, 0.000, b]]

    molecule = Structure(coordinates=coordinates,
                         symbols=symbols,
                         charge=chg,
                         multiplicity=mul)

    #- - -Perform RASCI QCHEM calculation and get information - - -
    print('\n# # # # # CALCULATIONS AT %.3f ANGS OF SEPARATION # # # # #' %b)
    rascidata, ee = ClcFs.RASCICalculation(molecule, bas, act, elec, occ,
                                         spin_mult, roots, verbose)
    #print(rascidata)
    print('\nInitial Active Space')
    print('Active electrons : ', elec)
    print('Active orbitals  : ', act)
    print('Occ orbitals     : ', occ)
    #- - - Extract information for calculation of natural orbitals - - -
    ovlp, mos, dm = wrp.rasci_ee_print(ee, False, False, False)

    #- - - Compute Natural orbitals and their occupations - - -
    print('- Python build-in NO coefficient matrix')
    nos, occs = oei.Build_NOs(dm, mos, ovlp)
    nos['qchem_order'] = ee['coefficients']['qchem_order']

    #- - - Redefine RAS2 based on occupations? - - -
#    occ_lev = 0.01
#    act_new, elec_new, occ_new, occ_larger = oei.NewRAS2(occs, occ_lev)
#    print('\nRedefined Active Space')
#    print('New active electrons : ', elec_new)
#    print('New active orbitals  : ', act_new)
#    print('New occ orbitals     : ', occ_new)
#    print('Lowest NO occ in RAS2: ', occ_larger)

    #- - - Avoid redefinition of RAS2 but cut orbital space there
    # Get first NO occupation outside RAS2
    occ_ras2, elec_ras2, occ_lower = oei.RAS2_occupations(occs, occ, elec, act)
    print('Occupations :', occs)
    print('occ_ras2 = ', occ_ras2)
    print('elec_ras2 = ', elec_ras2)
    print('occ_low = ', occ_lower)
    #- - - Compute the E>nu correlation energy with redefined active space and freezing virtuals
    # Calculation with redefined NewRAS2
    #dataFVirt, Elarger = ClcFs.RASFrozVirt(molecule, occs , occ_larger, nos, bas,
    #                     act_new, elec_new, occ_new, spin_mult, roots, verbose)

    #- - - Compute the E>nu correlation energy with redefined active space without freezing virtuals
    #dataFVirt, Elarger = ClcFs.RASFrozVirt(molecule, occs , 0.0, nos, bas,
    #            act_new, elec_new, occ_new, spin_mult, roots, verbose)


    #- - - Compute the E>nu correlation energy without changes in the RAS2 space, without freezing virtuals.
    dataFVirt, Elarger = ClcFs.RASFrozVirt(molecule, occs , 0.0, nos, bas,
                                           act, elec, occ, spin_mult, roots, verbose)


    #- - - Compute the E>nu correlation energy without changes in the RAS2 space, for first occupation
    # outside this space.
#    dataFVirt, Elarger = ClcFs.RASFrozVirt(molecule, occs , occ_lower, nos, bas,
#                                           act, elec, occ, spin_mult, roots, verbose)

    print('\nEnergy larger than threshold : ', Elarger[0])

    #- - - Compare with RASCI(h,p) total energy - - -
    pars_rasci = parser_rasci(rascidata)
    rasci_tot_e = pars_rasci['excited_states'][0]['total_energy']
    print('Total RASCI(h,p) energy : ', rasci_tot_e)

    print('\n# # # # # END OF CALCULATION AT %.3f ANGS OF SEPARATION # # # # #' %b)
    E_RAS.append(rasci_tot_e)
    E_RAS_NOs.append(Elarger[0][0])

# transform energy list in array
E_RAS = np.array(E_RAS)
E_RAS_NOs = np.array(E_RAS_NOs)
print(E_RAS)
print(E_RAS_NOs)
# print energies
print('\n# distance  RASCI(6,6)')
print('----------------------------')
for d, e in zip(distances, E_RAS):
    print('{:4.1f} : {:10.5f} '.format(d, e))

print('\n# distance  E_RAS_NOs_large')
print('----------------------------')
for d, e in zip(distances, E_RAS_NOs):
    print('{:4.1f} : {:10.5f} '.format(d, e))
