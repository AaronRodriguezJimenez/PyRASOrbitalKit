import numpy as np

# Here are functions to print and process the information from
# pyqchem information.

def rasci_ee_print(ee, OVLP, MO, RHO):
    """
    :param ee: electronic structure data disctionary from pyqchem
    param OVLP, MO, RHO are TRUE or FALSE, depending on what is requested
    :return: printing information and requested entity as a numpy array
    """
    ovlp = np.array(ee['overlap'])
    mos = np.array(ee['coefficients']['alpha'])
    tot_scf_dens = np.array(ee['total_scf_density'])
    #- - - Overlap matrix
    if OVLP == True:
        print('\n* * * * * * INFORMATION FROM NORMAL RASCI CALCULATION * * * * * * ')
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

    return ovlp, mos, tot_scf_dens


def AS_print(elec, act, occ):
    """
    :param elec: electrons in active space
    :param act:  active orbitals
    :param occ:  doubly occupied RAS1 orbitals
    :return:
    """
    print('\n- - - - - -  Active space - - - - - - ')
    print('-                                   - ')
    print('- Active electrons: %.d             - ' %elec)
    print('- Active orbitals: %.d              - ' %act)
    print('- Orbitals in RAS1: %.d             - ' %occ)
    print('-                                   - ')
    print('- - - - - - - - - - - - - - - - - - - ')

def NOONs_to_index(occs, threshold):
    """
    This function returns the number of large and small NOONs
    according to a threshold.
    :param occs: list of occupations
    :param threshold: numerical threshold
    :return: large_occs, small_occs
    """
    total_occs =  len(occs)
    for i, nu in enumerate(reversed(occs)):
        if nu < threshold:
            total_occs = total_occs - 1
        else:
            break

    large_occs = total_occs
    small_occs = len(occs) - total_occs

    return large_occs, small_occs

def NOONs_to_index_ver3(occs, threshold):
    """
    This function returns the number of large and small NOONs
    according to a threshold.
    :param occs: list of occupations
    :param threshold: numerical threshold
    :return: large_occs, small_occs whose sum is smaller than threshold
    """
    total_occs =  len(occs)
    sum1 = 0.00
    #sum2 = 0.00
    for i, nu in enumerate(reversed(occs)):
        sum1 = sum1 + nu
        if sum1 <= threshold:
            #sum2 = sum1
            total_occs = total_occs - 1
        else:
            break

    large_occs = total_occs
    small_occs = len(occs) - total_occs

    return large_occs, small_occs

def index_to_NOONs(occs, index, threshold):
    """
     This function returns the sum of occupations from a given set or occupations
    """
    total_occs = 0
    for i, nu in enumerate(reversed(occs)):
        if i < index:
            total_occs = total_occs + nu
        if nu < threshold:
            total_occs = total_occs + nu
        else:
            break

    return total_occs

def au2kcal(au):
    """
     This function converts Hartree to kcal/mol
    :param au: some input energy in hartree
    :return: energy in kcal/mol
    """
    kcal = 627.509
    return au*kcal

def au2ev(au):
    """
     This function converts Hartree to electron volts
    :param au: some input energy in hartree
    :return: energy in kcal/mol
    """
    ev = 27.211
    return au*ev