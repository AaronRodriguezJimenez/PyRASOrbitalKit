import math
from pyqchem.structure import Structure
import numpy as np
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_rasci import parser_rasci

#* * * * * * Functions * * * * * *
def OccsLOW(occslist, threshold):
    '''
    Function that recieves a list with al
    eliminates those larger than given th
    with those lower than the threshold.
    '''
    tail=[]
    num = 0 #Natural orbital number
    for nu in occslist:
        num = num+1
        if(nu < threshold):
            tail.append(nu)

    occ_num = len(tail)
    return tail, occ_num

def rasci_ee_print(ee, OVLP, MO, RHO):
    """
    :param ee: electronic structure data disctionary from pyqchem
    param OVLP, MO, RHO are TRUE or FALSE, depending on what is requested
    :return: printing information and requested entity as a numpy array
    """
    print('\n* * * * * * INFORMATION FROM NORMAL RASCI CALCULATION * * * * * * ')
    ovlp = np.array(ee['overlap'])
    mos = np.array(ee['coefficients']['alpha'])
    tot_scf_dens = np.array(ee['total_scf_density'])
    #- - - Overlap matrix
    if OVLP == True:
        print('- Overlap Matrix')
        print(ovlp)
    #- - - Alpha MO coefficient matrix
    if MO == True:
        print('- Alpha MO coefficient matrix')
        print(mos)
    #- - - Total SCF density
    if RHO == True:
        print('- Total SCF Density')
        print(tot_scf_dens)
    print('\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *')

    return ovlp, mos, tot_scf_dens

def Build_NOs(tot_scf_dens, mos, ovlp):
    """
    input:
    tot_scf_dens = total scf density matrix from electronic structure
    mos = molecular orbital coeficients from electronic structure
    ovlp = overlap matrix from electronic structure
    returns:
    natural orbital coefficient dictionary and numpy array of occupations
    """
    print('type of tot scf dens', type(tot_scf_dens))
    print('type mos', type(mos))
    print('type ovlp', type(ovlp))
    computed_NOs = {'alpha': [], 'beta': []}
    #- - - 1-PDM in MO basis
    print('\nTransform 1-PDM to AO basis')
    dens_mo = mos @ ovlp @ tot_scf_dens @ ovlp @ mos.T
    #- - - Diagonalization, get NOs and Occs
    weig, veig = np.linalg.eigh(dens_mo)
    # Sort eigenvalues and eigenvectors in descending order
    sorted_indices = np.argsort(weig)[::-1]
    weig = weig[sorted_indices]
    veig = veig[:, sorted_indices]
    # Compute natural orbitals in atomic orbital basis
    occupation_numbers = weig
    natural_orbitals_ao = np.array(veig.T @ mos, dtype=float)
    newnos = []
    for i in range(0, len(natural_orbitals_ao)):
        newnos.append(natural_orbitals_ao[i].tolist())

    print('\n- - - Creating ascending-order list of natural orbitals - - -')
    computed_NOs['alpha'] = newnos
    computed_NOs['beta'] = newnos
    print('type computed_NOs', type(computed_NOs))
    return computed_NOs, occupation_numbers

def RAS2_occupations(occs, occ, elec, act):
    """
    # Function that gives the NO occupations in RAS2
     Input: occs :  is a lists with the NO occupancies
            occ: RAS1 number of orbitals
            elec, act: RAS2 electrons and orbitals
     Output: occ_ras2 : NO occupations in RAS2
             elec_ras2 : electrons in active space
             occ_low: lower NO occupation in RAS2 equivalent
                      to occ_larger in NewRAS2 function
    """
    print("\nGetting active space NO occupations")
    print('original electrons in ras2 =', elec)
    id_occ_ras1 = occ
    id_occ_ras2 = occ + act
    lst_occ_ras2 = np.arange(id_occ_ras1,id_occ_ras2)
    print('lst_occ_ras2 : ', lst_occ_ras2)
    occ_ras2 = []
    for i in lst_occ_ras2:
        print('id_occs_in_ras2', i , occs[i])
        occ_ras2.append(occs[i])

    elec_ras2 = math.ceil(sum(occ_ras2))
    occ_low = occ_ras2[len(occ_ras2)-1]
    print('determined electrons in ras2 =', elec_ras2)
    return occ_ras2, elec_ras2, occ_low

def RASCICalculation(molecule, bas, act, elec, occ, spin_mult, roots, verbose):
    """
     This function performs a standard RASCI(h,p) calculation:
     molecule = pyqchem object with coordinates and molecule info
     bas = choosen basis set
     act = RAS2 active orbitals
     elec = RAS2 active electrons
     occ = RAS1 ocupied orbitals
     spin_mult = RAS roots multiplicity
     roots = number of roots to compute
     verbose = RAS_PRINT keyword value, printing level

     Returns:
     data1  = output from QCHEM
     ee = electronic structure information
    """

    #- - - RASCI Normal Calculation - - -
    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          exchange='hf',
                          unrestricted=False,
                          correlation='rasci',
                          basis=bas,
                          ras_act=act,
                          ras_elec=elec,
                          ras_occ=occ,
                          ras_spin_mult=spin_mult,
                          ras_roots=roots,
                          ras_do_part=True,
                          ras_do_hole=True,
                          max_cis_cycles=999,
                          set_iter=999,
                          mem_total=2000,
                          mem_static=100,
                          symmetry=False,
                          ras_print=verbose,
                          state_analysis=True)


    data1, ee = get_output_from_qchem(qc_input,
                                      processors=4,
                                      force_recalculation=False,
                                      return_electronic_structure=True,
                                      store_full_output=True)
    return data1, ee
def RASFrozVirt(molecule, occs, thrs, NOsguess,
                bas, act, elec, occ, spin_mult, roots, verbose):
    """
     This function performs a RASCI calculation for the E>_nu with
     computed NOs as initial guess
     returns:
     data  = output from QCHEM
     eNOF = DFT energy list computed for each root
    """
    #- - - First determine the virtuals to be frozen
    occs_lower , frozenorb = OccsLOW(occs, thrs)
    print('Occupations lower than occupation %.4f threshold :' %(thrs),
          occs_lower)
    print('Frozen virtual orbitals:', frozenorb)
    #- - - RASCI frozen virtual Calculation - - -
    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          exchange='hf',
                          unrestricted=False,
                          correlation='rasci',
                          basis=bas,
                          scf_guess=NOsguess,
                          max_scf_cycles=0,
                          ras_act=act,
                          ras_elec=elec,
                          ras_occ=occ,
                          ras_spin_mult=spin_mult,
                          ras_roots=roots,
                          set_iter=999,
                          ras_do_part=True,
                          ras_do_hole=True,
                          max_cis_cycles=999,
                          mem_total=2000,
                          mem_static=100,
                          symmetry=False,
                          n_frozen_virt=frozenorb,
                          ras_print=verbose,
                          extra_rem_keywords={'NEW_AUX_ORDER':1})

    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 force_recalculation=True,
                                 return_electronic_structure=False,
                                 store_full_output=True)

    parsed_data = parser_rasci(data)

    energies = []
    energies.append([state['total_energy'] for state in parsed_data['excited_states']])

    return data, energies

#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#- - - QCHEM KEYWORD VALUES FOR RASCI CALCULATION - - -
bas = '6-31g*'
act = 6
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
               [0.000, 0.000, 1.00]]

molecule = Structure(coordinates=coordinates,
                     symbols=symbols,
                     charge=chg,
                     multiplicity=mul)

#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
E_RAS = []
E_RAS_NOs = []
#- - -Perform RASCI QCHEM calculation and get information - - -
rascidata, ee = RASCICalculation(molecule, bas, act, elec, occ,
                                       spin_mult, roots, verbose)

print(rascidata)
#- - - Extract information for calculation of natural orbitals - - -
ovlp, mos, dm = rasci_ee_print(ee, False, False, False)

#- - - Compute Natural orbitals and their occupations - - -
ovlp = np.array(ovlp)
nos, occs = Build_NOs(dm, mos, ovlp)

print('NOs: ')
print(nos)
print('NOONs: ')
print(occs)

import pyrasorbitalkit.OrbFs as OrbFs
#Compare with other NOs:
nos, occs = OrbFs.Build_NOs(dm, mos, ovlp)
print('NOs 2: ')
print(nos)
print('NOONs 2: ')
print(occs)



#- - - Avoid redefinition of RAS2 but cut orbital space there
# Get first NO occupation outside RAS2
occ_ras2, elec_ras2, occ_lower = RAS2_occupations(occs, occ, elec, act)

# Calculation keeping RAS2 fixed with frozen virtuals equal to occ_lower
dataFVirt, Elarger = RASFrozVirt(molecule, occs , occ_lower, nos, bas,
                                  act, elec, occ, spin_mult, roots, verbose)

#- - - Compare with RASCI(h,p) total energy - - -
pars_rasci = parser_rasci(rascidata)
rasci_tot_e = pars_rasci['excited_states'][0]['total_energy']
E_RAS.append(rasci_tot_e)
E_RAS_NOs.append(Elarger[0][0])

#- - - Compare with RASCI(h,p) total energy - - -
pars_rasci = parser_rasci(rascidata)
rasci_tot_e = pars_rasci['excited_states'][0]['total_energy']

# transform energy list in array
E_RAS = np.array(E_RAS)
E_RAS_NOs = np.array(E_RAS_NOs)

# print energies
for e1, e2 in zip(E_RAS, E_RAS_NOs):
    print('\nRASCI(6,6)   RAS_with_NOs_guess')
    print('{:10.5f} : {:10.5f} '.format(e1, e2))
