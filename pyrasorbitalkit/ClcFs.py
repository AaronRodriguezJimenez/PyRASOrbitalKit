import numpy as np
from pyqchem.qchem_core import get_output_from_qchem
from pyqchem.qc_input import QchemInput
from pyqchem.parsers.parser_rasci import parser_rasci
#from parser_rasci_rasnof0 import parser_rasci
from pyrasorbitalkit import OrbFs
from pyrasorbitalkit import wrp


def RHF_calc(molecule, bas):
    # - - - - RHF calculation - - - -
    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          method='hf',
                          unrestricted=False,
                          basis=bas,
                          set_iter=999,
                          mem_total=2000,
                          mem_static=100,
                          symmetry=False,
                          sym_ignore=False)

    data, ee = get_output_from_qchem(qc_input,
                                     processors=4,
                                     force_recalculation=False,
                                     return_electronic_structure=True,
                                     store_full_output=True)
    return data, ee

def RHF_calc_guess_bas2(molecule, bas, bastwo, mos):
    # - - - - RHF calculation - - - -
    # Useful for concatenated calcualtions within QCHEM basis set projection
    #
    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          method='hf',
                          unrestricted=False,
                          scf_guess=mos,
                          basis=bas,
                          basis2=bastwo,
                          max_scf_cycles=999,
                          set_iter=999,
                          mem_total=2000,
                          mem_static=100,
                          symmetry=False,
                          sym_ignore=False)

    data, ee = get_output_from_qchem(qc_input,
                                     processors=4,
                                     force_recalculation=False,
                                     return_electronic_structure=True,
                                     store_full_output=True)
    return data, ee


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

    # - - - RASCI Normal Calculation - - -
    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          exchange='hf',
                          unrestricted=False,
                          correlation='rasci',
                          basis=bas,
                          ras_act=act,
                          ras_elec=elec,
                          max_scf_cycles=999,
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
                          sym_ignore=False,
                          ras_print=verbose,
                          state_analysis=True)

    data1, ee = get_output_from_qchem(qc_input,
                                      processors=4,
                                      force_recalculation=False,
                                      return_electronic_structure=True,
                                      store_full_output=True)
    return data1, ee


def RASCIGuess(molecule, mos, bas, act, elec, occ, hole, part, spin_mult, roots, verbose):
    """
     This function performs a standard RASCI(h,p) calculation useful for
     concatenated calculations:
     molecule = pyqchem object with coordinates and molecule info
     mos = guess orbitals
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

    # - - - RASCI Normal Calculation - - -
    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          exchange='hf',
                          unrestricted=False,
                          correlation='rasci',
                          scf_guess=mos,
                          basis=bas,
                          ras_act=act,
                          ras_elec=elec,
                          max_scf_cycles=999,
                          ras_occ=occ,
                          ras_spin_mult=spin_mult,
                          ras_roots=roots,
                          ras_do_part=part,
                          ras_do_hole=hole,
                          max_cis_cycles=999,
                          set_iter=999,
                          mem_total=2000,
                          mem_static=100,
                          symmetry=False,
                          sym_ignore=False,
                          ras_print=verbose,
                          state_analysis=True)

    data1, ee = get_output_from_qchem(qc_input,
                                      processors=4,
                                      force_recalculation=False,
                                      return_electronic_structure=True,
                                      store_full_output=True)
    return data1, ee


def RASCIGuess_bas2(molecule, mos, bas, bastwo, act, elec, occ, hole, part, spin_mult, roots, verbose):
    """
     This function performs a standard RASCI(h,p) calculation useful for
     concatenated calculations, uses BASIS2 keyword and guess from aprevious
     calcultion within a smaller basis set:
     molecule = pyqchem object with coordinates and molecule info
     mos = guess orbitals
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

    # - - - RASCI Normal Calculation - - -
    qc_input = QchemInput(molecule,
                          jobtype='sp',
                          exchange='hf',
                          unrestricted=False,
                          correlation='rasci',
                          scf_guess=mos,
                          basis=bas,
                          basis2=bastwo,
                          ras_act=act,
                          ras_elec=elec,
                          max_scf_cycles=999,
                          ras_occ=occ,
                          ras_spin_mult=spin_mult,
                          ras_roots=roots,
                          ras_do_part=part,
                          ras_do_hole=hole,
                          max_cis_cycles=999,
                          set_iter=999,
                          mem_total=2000,
                          mem_static=100,
                          symmetry=False,
                          sym_ignore=False,
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
     occs: list of NOONs in numpy format
     thrs: Threshold for occupations
     returns:
     data  = output from QCHEM
     eNOF = DFT energy list computed for each root
    """
    # - - - For the case of a high-spin ref state, add charge to avoid error
    # - - - this is not supposed to affect the results.
    if np.mod(molecule.multiplicity, 2) == 0:
        molecule.charge = -1

    molecule.multiplicity = 1

    # - - - First determine the virtuals to be frozen
    occs_lower, frozenorb = OrbFs.OccsLOW(occs, thrs)
    # print('Occupations lower than occupation %.4f threshold :' %(thrs),
    #      occs_lower)
    print('Frozen virtual orbitals:', frozenorb)
    if (len(occs) - frozenorb) < occ + act:
        print(len(occs) - frozenorb, '//', occ + act)
        print('\n***** ERROR: THRESHOLD TOO LARGE, E>nu CAN NOT BE COMPUTED *******')
        print('Change occupation threshold, here is the list of occupations')
        print(occs)
        print('Returning values : data = 0, ee = 0, energies = 0')
        return 0, 0, 0
    else:
        # - - - RASCI frozen virtual Calculation - - -
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
                              sym_ignore=False,
                              n_frozen_virt=frozenorb,
                              ras_print=verbose,
                              extra_rem_keywords={'NEW_AUX_ORDER': 1})

        data, ee = get_output_from_qchem(qc_input,
                                         processors=4,
                                         force_recalculation=True,
                                         return_electronic_structure=True,
                                         store_full_output=True)

        parsed_data = parser_rasci(data)

        energies = []
        energies.append([state['total_energy'] for state in parsed_data['excited_states']])

        return data, ee, energies


def RASFrozVirt3(molecule, occs, thrs, NOsguess,
                 bas, act, elec, occ, spin_mult, roots, verbose):
    """
     This function performs a RASCI calculation for the E>_nu with
     computed NOs as initial guess (For RAS-NOF3 version only)
     returns:
     data  = output from QCHEM
     eNOF = DFT energy list computed for each root
    """
    # - - - For the case of a high-spin ref state, add charge to avoid error
    # - - - this is not supposed to affect the results.
    if np.mod(molecule.multiplicity, 2) == 0:
        molecule.charge = -1

    molecule.multiplicity = 1

    # - - - First determine the virtuals to be frozen
    occs_lower, frozenorb = OrbFs.SumOccsLow(occs, thrs)
    print('Occupations lower than occupation %.4f threshold :' % (thrs),
          occs_lower)
    print('Frozen virtual orbitals:', frozenorb)
    if (len(occs) - frozenorb) < occ + act:
        print(len(occs) - frozenorb, '//', occ + act)
        print('\n***** ERROR: THRESHOLD TOO LARGE, E>nu CAN NOT BE COMPUTED *******')
        print('Change occupation threshold, here is the list of occupations')
        print(occs)
        print('Returning values : data = 0, ee = 0, energies = 0')
        return 0, 0, 0
    else:
        # - - - RASCI frozen virtual Calculation - - -
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
                              sym_ignore=False,
                              n_frozen_virt=frozenorb,
                              ras_print=verbose,
                              extra_rem_keywords={'NEW_AUX_ORDER': 1})

        data, ee = get_output_from_qchem(qc_input,
                                         processors=4,
                                         force_recalculation=True,
                                         return_electronic_structure=True,
                                         store_full_output=True)

        parsed_data = parser_rasci(data)

        energies = []
        energies.append([state['total_energy'] for state in parsed_data['excited_states']])

        return data, ee, energies


def RASNOF0(molecule, occthrs, version, nofdft, bas, act, elec, occ,
            spin_mult, roots, verbose):
    """
     This function performs a RAS-NOF1 calculation for the E_dft:
     Returns:
     eNOF = DFT energy list computed for each root
    """

    # - - - RAS-NOF on-top Calculation - - -
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
                          sym_ignore=False,
                          ras_natorb=True,
                          ras_print=verbose,
                          extra_rem_keywords={'ras_nof': True,
                                              'ras_nof_thresh_occ': occthrs,
                                              'ras_nof_dft ': nofdft,
                                              'ras_nof_version': version})

    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 force_recalculation=True,
                                 return_electronic_structure=False,
                                 store_full_output=True)

    parsed_data = parser_rasci(data)

    energies = []
    nof_energies = []
    energies.append([state['total_energy'] for state in parsed_data['excited_states']])
    nof_energies.append([state['total_nof_energy'] for state in parsed_data['excited_states']])
    return data, energies, nof_energies


def RASNOF30_sum(molecule, occthrs, version, nofdft, bas, act, elec, occ,
                 spin_mult, roots, verbose):
    """
     This function performs a RAS-NOF1 calculation for the E_dft:
     Returns:
     eNOF = DFT energy list computed for each root
    """

    # - - - RAS-NOF on-top Calculation - - -
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
                          sym_ignore=False,
                          ras_print=verbose,
                          extra_rem_keywords={'ras_nof': True,
                                              'ras_nof_thresh_sum': occthrs,
                                              'ras_nof_dft ': nofdft,
                                              'ras_nof_version': version})

    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 force_recalculation=True,
                                 return_electronic_structure=False,
                                 store_full_output=True)

    parsed_data = parser_rasci(data)

    energies = []
    nof_energies = []
    energies.append([state['total_energy'] for state in parsed_data['excited_states']])
    nof_energies.append([state['total_nof_energy'] for state in parsed_data['excited_states']])
    return data, energies, nof_energies


def NOFDFT1(molecule, occthrs, bas, act, elec, occ, spin_mult, roots, verbose):
    """
     This function performs a RAS-NOF1 calculation for the E_dft:
     Returns:
     eNOF = DFT energy list computed for each root
    """

    # - - - RAS-NOF1 (Savin) Calculation - - -
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
                          sym_ignore=False,
                          ras_print=verbose,
                          extra_rem_keywords={'ras_nof': True,
                                              'ras_nof_thresh_occ': occthrs,
                                              'ras_nof_dft ': 'VWN1_NOF',
                                              'ras_nof_version': 1})

    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 force_recalculation=True,
                                 return_electronic_structure=False,
                                 store_full_output=True)

    parsed_data = parser_rasci(data)

    energies = []
    nof_energies = []
    energies.append([state['total_energy'] for state in parsed_data['excited_states']])
    nof_energies.append([state['total_nof_energy'] for state in parsed_data['excited_states']])
    eNOF = []
    print(nof_energies[0][0] - energies[0][0])
    for e, nof in zip(energies, nof_energies):
        eNOF.append(nof[0] - e[0])  # Only works for one root NEED furhter revision

    return eNOF


def NOFDFT2_savin(molecule, occthrs, bas, act, elec, occ, spin_mult, roots, verbose):
    """
     This function performs a RAS-NOF1 calculation for the E_dft:
     Returns:
     eNOF = DFT energy list computed for each root
    """

    # - - - RAS-NOF12 (Savin) Calculation - - -
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
                          sym_ignore=False,
                          ras_print=verbose,
                          extra_rem_keywords={'ras_nof': True,
                                              'ras_nof_thresh_occ': occthrs,
                                              'ras_nof_dft ': 'VWN1_NOF2',
                                              'ras_nof_version': 2})

    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 force_recalculation=True,
                                 return_electronic_structure=False,
                                 store_full_output=True)

    parsed_data = parser_rasci(data)

    energies = []
    nof_energies = []
    energies.append([state['total_energy'] for state in parsed_data['excited_states']])
    nof_energies.append([state['total_nof_energy'] for state in parsed_data['excited_states']])
    eNOF = []
    print(nof_energies[0][0] - energies[0][0])
    for e, nof in zip(energies, nof_energies):
        eNOF.append(nof[0] - e[0])  # Only works for one root NEED furhter revision

    return eNOF


def NOFDFT2_gori(molecule, occthrs, bas, act, elec, occ, spin_mult, roots, verbose):
    """
     This function performs a RAS-NOF1 calculation for the E_dft:
     Returns:
     eNOF = DFT energy list computed for each root
    """

    # - - - RAS-NOF2 (Gori-Giorgi) Calculation - - -
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
                          sym_ignore=False,
                          ras_print=verbose,
                          extra_rem_keywords={'ras_nof': True,
                                              'ras_nof_thresh_occ': occthrs,
                                              'ras_nof_dft': 'VWN1_NOF3',
                                              'ras_nof_version': 2})

    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 force_recalculation=True,
                                 return_electronic_structure=False,
                                 store_full_output=True)

    parsed_data = parser_rasci(data)

    energies = []
    nof_energies = []
    energies.append([state['total_energy'] for state in parsed_data['excited_states']])
    nof_energies.append([state['total_nof_energy'] for state in parsed_data['excited_states']])
    eNOF = []
    print(nof_energies[0][0] - energies[0][0])
    for e, nof in zip(energies, nof_energies):
        eNOF.append(nof[0] - e[0])  # Only works for one root NEED furhter revision

    return eNOF


def NOFDFT3_savin(molecule, occs, sumthrs, bas, act, elec, occ, spin_mult, roots, verbose):
    """
     This function performs a RAS-NOF1 calculation for the E_dft:
     Returns:
     eNOF = DFT energy list computed for each root
    """
    thresh_sum_list, froz_virt = OrbFs.SumOccsLow(occs, sumthrs)
    print('- - - Orbitals lower than threshold in RAS-NOF3: ', thresh_sum_list)
    print('- - - These orbitals sum in total to: ', np.sum(thresh_sum_list))
    print('- - - which must be lower than the sum threshold: ', sumthrs)
    print('- - - total number or orbitals with occupations lower than sum threshold: ', froz_virt)
    # - - - RAS-NOF3 (Savin) Calculation - - -
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
                          sym_ignore=False,
                          ras_print=verbose,
                          extra_rem_keywords={'ras_nof': True,
                                              'ras_nof_thresh_sum': np.sum(thresh_sum_list),  # sumthrs,
                                              'ras_nof_dft': 'VWN1_NOF2',
                                              'ras_nof_version': 3})

    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 force_recalculation=True,
                                 return_electronic_structure=False,
                                 store_full_output=True)

    parsed_data = parser_rasci(data)

    energies = []
    nof_energies = []
    energies.append([state['total_energy'] for state in parsed_data['excited_states']])
    nof_energies.append([state['total_nof_energy'] for state in parsed_data['excited_states']])
    eNOF = []
    print(nof_energies[0][0] - energies[0][0])
    for e, nof in zip(energies, nof_energies):
        eNOF.append(nof[0] - e[0])  # Only works for one root NEED furhter revision

    return eNOF


def NOFDFT3_gori(molecule, occs, sumthrs, bas, act, elec, occ, spin_mult, roots, verbose):
    """
     This function performs a RAS-NOF1 calculation for the E_dft:
     Returns:
     eNOF = DFT energy list computed for each root
    """
    thresh_sum_list, froz_virt = OrbFs.SumOccsLow(occs, sumthrs)
    print('- - - Orbitals lower than threshold in RAS-NOF3: ', thresh_sum_list)
    print('- - - These orbitals sum in total to: ', np.sum(thresh_sum_list))
    print('- - - which must be lower than the sum threshold: ', sumthrs)
    print('- - - total number or orbitals with occupations lower than sum threshold: ', froz_virt)
    # - - - RAS-NOF3 (Gori-Giorgi) Calculation - - -
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
                          sym_ignore=False,
                          ras_print=verbose,
                          extra_rem_keywords={'ras_nof': True,
                                              'ras_nof_thresh_sum': np.sum(thresh_sum_list),  # sumthrs,
                                              'ras_nof_dft': 'VWN1_NOF3',
                                              'ras_nof_version': 3})

    data = get_output_from_qchem(qc_input,
                                 processors=4,
                                 force_recalculation=True,
                                 return_electronic_structure=False,
                                 store_full_output=True)

    parsed_data = parser_rasci(data)

    energies = []
    nof_energies = []
    energies.append([state['total_energy'] for state in parsed_data['excited_states']])
    nof_energies.append([state['total_nof_energy'] for state in parsed_data['excited_states']])
    eNOF = []
    print(nof_energies[0][0] - energies[0][0])
    for e, nof in zip(energies, nof_energies):
        eNOF.append(nof[0] - e[0])  # Only works for one root NEED furhter revision

    return eNOF


def RAS_NO(molecule, NOsguess, bas, act, elec, occ,
           spin_mult, roots, verbose):
    """
     This function performs a RASCI calculation for the E>_nu with
     computed NOs as initial guess
     NOsguess = Natural Orbital coefficients guess
     returns:
     data  = output from QCHEM
     ee =  electronic structure information from calculation
     energies = total energies from states in calculation
    """
    # - - - For the case of a high-spin ref state, add charge to avoid error
    # - - - this is not supposed to affect the results.
    if np.mod(molecule.multiplicity, 2) == 0:
        molecule.charge = -1

    molecule.multiplicity = 1

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
                          sym_ignore=False,
                          ras_print=verbose,
                          state_analysis=True,
                          extra_rem_keywords={'NEW_AUX_ORDER': 1})

    data, ee = get_output_from_qchem(qc_input,
                                     processors=4,
                                     force_recalculation=True,
                                     return_electronic_structure=True,
                                     store_full_output=True)

    parsed_data = parser_rasci(data)

    energies = []
    energies.append([state['total_energy'] for state in parsed_data['excited_states']])

    return data, ee, energies


def nof_sc(molecule, NOsguess, tolerance, bas, act, elec, occ,
           spin_mult, roots, verbose):
    """
    tolerance - is the tolerance for the SC procedure
    NOsguess - are the NO coefficient matrix guess

    This function computes the self-consistent energy using NOs as input
    :return: RAS(sc-nof) energy computed with the complete RAS3 space
    """

    # - - - Initial NOF calc for comparison, ref energ for the SC:
    data, eeNO, energies = RAS_NO(molecule, NOsguess, bas, act, elec, occ,
                                  spin_mult, roots, verbose)

    # - - - Get Natural orbitals and occupations
    print('- - -SC procedure - Electronic structure NOONs - - -')
    nos = eeNO['nato_coefficients']
    nos['qchem_order'] = eeNO['coefficients']['qchem_order']
    occs = eeNO['nato_occupancies']['alpha']
    if verbose > 2:
        print(occs)

    ras_no_e = energies[0][0]
    print('\n* Initial RAS(NO) energy : ', ras_no_e)

    # - - -Initialize the SC procedure:
    max_iterations = 100
    for iteration in range(max_iterations + 1):
        # Original variables from previous cycle, for comparison
        E_old = ras_no_e
        nos_old = nos

        # Perform RAS-NO calculation:
        # - - - (h,p) total scheme - - -
        data, eeNO, energies = RAS_NO(molecule, nos_old, bas, act, elec, occ,
                                      spin_mult, roots, verbose)
        Eel = energies[0][0]
        # - - - (h,n) frozen virtual scheme - - - (pending)
        # - - - (h) complete virtual space frozen - - - (pending)

        # Compute/get new natural orbitals
        #ovl, mos, dm = wrp.rasci_ee_print(eeNO, False, False, False)
        # nos, occs = OrbFs.Build_NOs(dm, mos, ovl)
        # nos['qchem_order'] = eeNO['coefficients']['qchem_order']
        print('Electronic structure NOONs')
        nos = eeNO['nato_coefficients']
        nos['qchem_order'] = eeNO['coefficients']['qchem_order']
        occs = eeNO['nato_occupancies']['alpha']

        print('- - - NOF-SC Procedure - - -')
        print(f"- Iteration : ", iteration + 1)
        print(f"- RAS(NO) Energy root 1: ", Eel)
        if (np.fabs(Eel - E_old) < tolerance) and (iteration > 0):
            print("- - - Convergence achieved - - -")
            occ, elec, occ_low = OrbFs.RAS2_occupations(occs, occ, elec, act)
            print("Lower occupation in RAS2 :", occ_low)
            print('Final Energy Calc with NOs guess: ', energies[0][0])
            print('NO Occs :', occs)
            print(data)
            break

        if iteration == max_iterations - 1:
            print("\nMaximum iterations reached, calculation did not converged.")

    return energies[0][0]