"""
Author: tzs820, s240670, Emil Henningsen
Module for constructing different Hamiltonians used for describing atomic arrays interacting with light.
"""

def coherence_operators(i, j, N):
    """
    Construct the space of coherence operators at i & j and identity elsewhere.
    Ref.: Asenjo-Garcia et al. equation (5)

    --- Parameters

    i: integer, i'th subspace of raising operator
    j: integer, j'th subspace of lowering operator
    N: integer, number of subspaces

    --- Return

    space: Qobj, array (2^N x 2^N) of 0's and 1's representing constructed space.
    """
    from qutip import Qobj, tensor, qeye

    sigma_ge = Qobj([[0, 1], [0, 0]])                         #deexcitation
    sigma_eg = Qobj([[0, 0], [1, 0]])                         #excitation

    if i == 0:
        space = sigma_eg
    elif j == 0:
        space = sigma_ge
    else:
        space = qeye(2)
    for k in range(1, N):
        if k == i:
            space = tensor([space, sigma_eg])
        elif k == j:
            space = tensor([space, sigma_ge])
        else:
            space = tensor([space, qeye(2)])
        #print(k)

    return space

def H_eff(N, w0, D, G):
    """
    Function for calculating effective Hamiltonian, ref.: Asenjo-Garcia et al. equation (5)

    --- Parameters:

    N: integer, number of dipoles (subspaces)
    w0: float, transition frequency of dipoles
    D: array (1 x 3) of complex floats, column-vector of dipole transition elements
    G: array (N x N x 3 x 3) of complex floats, Green's Tensors for a given lattice of dipoles.

    --- Return:

    H_eff: Qobj, array (2^N x 2^N) of complex floats
    """
    from qutip import Qobj
    from scipy.constants import mu_0

    H_eff = 0                                                       #Effective Hamiltonian of the system
    for i in range(N):                                              #Fill the effective Hamiltonian
        for j in range(N):
            if i == j:
                continue
            else:
                H_eff += (-mu_0 * w0**2) * D.trans() * Qobj(G[i,j]) * D * coherence_operators(i, j, N)
    return H_eff

def H(N, w0, H_eff):
    """
    Function for calculating full Hamiltonian, ref.: Asenjo-Garcia et al. equation (5)

    --- Parameters:

    N: integer, number of dipoles (subspaces)
    w0: float, transition frequency of dipoles
    H_eff: array (2^N x 2^N) of complex floats, effective Hamiltonian

    --- Return:

    H: Qobj, array (2^N x 2^N) of complex floats, full Hamiltonian.
    """
    from qutip import Qobj, qeye, tensor
    from scipy.constants import hbar

    sigma_ee = Qobj([[0, 0], [0, 1]])                             #excited state
    H = 0
    for i in range(N):

        if i == 0:
            space = sigma_ee
        else:
            space = qeye(2)
        for k in range(1, N):
            if k == i:
                space = tensor([space, sigma_ee])
            else:
                space = tensor([space, qeye(2)])

        H += hbar * w0 * space
        
    H += H_eff
    return H

def block(N, G, n):
    """
    Desciption: TODO
    """
    from numpy import zeros

    #If only one directional vector is passed, e.g. ez, then fill (N x 3) dim array of given direction
    if n.ndim == 1:
        arr = zeros((N,3))
        for i in range(N):
            arr[i] = n
        n = arr
    
    #We want Hamiltonian in basis of {|e_j>}, meaning one excitation on j'th subspace.
    #Do this by computing every matrix element of NxN matrix. It is the block of one excitation in full Hamiltonian
    block = zeros((N,N), dtype=complex)

    for i in range(N):
        for j in range(N):
            if i == j:
                block[i, j] += n[i].transpose() @ G[i,j] @ n[i] #WARNING: This might be wrong?!?
            else:
                block[i, j] += n[i].transpose() @ G[i,j] @ n[i]        #Dipoles polarized along z-direction.  

    return block

def scalarham(N, rij, w0=1, dimensionless=True):
    """
    See scalar() under GreensTensor.py
    """
    from GreensTensor import scalar

    h = scalar(N, rij, w0, dimensionless)
    return h