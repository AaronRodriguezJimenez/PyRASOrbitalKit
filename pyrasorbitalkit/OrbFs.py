import numpy as np
import scipy
from functools import reduce

#- - - FUNCTIONS FOR handling orbitals and active space - - -
#
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


def SumOccsLow(occslist, sumthresh):
    '''
    Function for RAS-NOF3, computes the s
    threshold for the sum, and gives the
    to be frozen.
    '''
    #print(type(occslist))
    nt = len(occslist)
    for i in range(nt-1):
        sp = np.sum(occslist[-i-1:]) #sum
        if sp > sumthresh:
            return occslist[nt-i:], i
    return [] , 0


def MAxValList(lst):
    max_val = max(lst)
    max_idx = lst.index(max_val)
    return max_idx, max_val


def ReorderNOs(mos, nos, ovlp):
    """
    Order based on Maximum overlap with MOs
    mosT = MO coefficients transposed vector
    max_j_id = the label for the NO with maximum OVLP with
               a given MO
    max_UdagSU = maximum value of the OVLP between MO and NO
    USN_list = list of ovelaps of a given MO with all NOs
    output: reordered dictionary of NOs based on maximun ovlp with MOs
    """
    newnos = []
    ordered_NOs = {'alpha': [], 'beta': []}
    mos = np.array(mos)
    nos = np.array(nos)
    for i in range(0, len(mos)):
        USN_list = []
        for j in range(0, len(nos)):
                mosT = mos[i].reshape(1,-1)
                UdagSU = mosT @ ovlp @ nos[j]
                USN_list.append(abs(UdagSU))

        max_j_id, max_UdagSU = MAxValList(USN_list)
#        print('i,  max_j_id,  max_UdagSU')
#        print(i, max_j_id, max_UdagSU)
        newnos.append(nos[max_j_id].tolist())
    print('\n- - - Creating ordered list of natural orbitals - - -')
    ordered_NOs['alpha'] = newnos
    ordered_NOs['beta'] = newnos

    return ordered_NOs

def NewRAS2(occs, level, max_tries=10):
    """
    #- - - Following J. Mol. Struct. (Theochem) 434 (1998) 239-245
     Input: occs :  is a lists with the NO occupancies
            level : is a threshold value to discriminate occupations
     Output: act_new : new active orbitals in RAS2
             elec_new : electrons in active space
             occ_new : mew doubly occupied orbitals in RAS1
             occ_larger: NO occupation at the top of the new RAS2
    If there are zero electrons or zero orbitals in RAS2 the selection
    bounds are readjusted.
    """
    tolerance = 0.000000001
    upper_val = 2.00 - level
    lower_val = level
    print("\nDefining second new active space based on NO occupations")
    #print('occupations :', occs)
    print("\nRAS2_OCC defined for NOs with occupations between %.5f-%.5f"
          % (upper_val, lower_val) )

    active_occs = []
    ras1_elec = 0
    nos_active_id = []

    for i, occ in enumerate(occs):
        if occ <= (2-level) and occ >= level:
           #print('occ in RAS2: ', occ)
           active_occs.append(occ)
           nos_active_id.append(i)
           #print('nos_active_id: ' ,nos_active_id)
           #print('len_nos_active_id: ', len(nos_active_id))
        if occ > (2-level):
           ras1_elec += 2

    act_new = len(active_occs)
    elec_new = sum(active_occs)
    elec_new = int(round(elec_new))
    occ_new = int(round(ras1_elec)/2)
   #if act_new == 0 or elec_new == 0:
    if not nos_active_id or elec_new ==0:
       print('\nNew RAS2 contains zero orbitals')
       for i in range(max_tries):
           print('Reducing occ_lev by half its original value')
           level = level-0.0001
           if level >= tolerance:
               act_new, elec_new, occ_new, occ_larger = NewRAS2(occs, level, max_tries)
           if act_new > 0 and elec_new > 0:
               return act_new, elec_new, occ_new, occ_larger

    id_larger = occ_new + act_new -1
    occ_larger = occs[id_larger]
    #print('Larger NO occupation: ' ,occ_larger)

    return act_new, elec_new, occ_new, occ_larger

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
    print('Original electrons in RAS2 from input =', elec)
    id_occ_ras1 = occ
    id_occ_ras2 = occ + act
    lst_occ_ras2 = np.arange(id_occ_ras1,id_occ_ras2)
    print('lst_occ_ras2 : ', lst_occ_ras2)
    occ_ras2 = []
    for i in lst_occ_ras2:
        print('id_occs_in_ras2', i , occs[i])
        occ_ras2.append(occs[i])

    elec_ras2 = sum(occ_ras2)
    occ_low = occ_ras2[len(occ_ras2)-1]
    print('Determined electrons in RAS2 from NOONs =', elec_ras2)
    return occ, elec, occ_low


def Build_NOs(tot_scf_dens, mos, ovlp):
    """
    input:
    tot_scf_dens = total scf density matrix from electronic structure
    mos = molecular orbital coeficients from electronic structure
    ovlp = overlap matrix from electronic structure
    returns:
    natural orbital coefficient dictionary and numpy array of occupations
    """
    print('- - - BUILDING NOs AND NOONs with Python - - -')
    computed_NOs = {'alpha': [], 'beta': []}
    # - - - 1-PDM in MO basis
    #print('\nTransform 1-PDM to AO basis')
    dens_mo = mos @ ovlp @ tot_scf_dens @ ovlp @ mos.T
    # - - - Diagonalization, get NOs and Occs
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

    #print('\n- - - Creating ascending-order list of natural orbitals - - -')
    computed_NOs['alpha'] = newnos
    computed_NOs['beta'] = newnos
    return computed_NOs, occupation_numbers

    # print('\n- - - Creating ascending-order list of natural orbitals - - -')
    computed_NOs['alpha'] = newnos
    computed_NOs['beta'] = newnos
    # print('type computed_NOs', type(computed_NOs))
    return computed_NOs, occupation_numbers