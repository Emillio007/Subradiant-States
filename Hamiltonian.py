"""
Author: tzs820, s240670, Emil Henningsen
Module for constructing different Hamiltonians used for describing atomic arrays interacting with light.
"""

from numpy import ndarray

class Hamiltonian:

    #flags:
    initialized = False
    decomposed = False

    N = None
    hamiltonian = None
    eigvec = None
    eigval = None

    #Standard exception:
    e = Exception("Hamiltonian is empty. ")

    def __init__(self, N : int = None, ham : ndarray = None):
        if ham == None:
            self.initialized = False
        else:
            self.initialized = True
        self.hamiltonian = ham
        self.N = N

    #Assertations:
    def assertInit(self):
        if self.initialized == False:
            raise self.e
        else:
            pass

    """
    GET methods:
    """

    def getN(self) -> int:
        return self.N

    def getHam(self) -> ndarray:
        return self.hamiltonian
    

    #Specialized:

    def getEigenDecomp(self, sort : bool = True) -> tuple[ndarray, ndarray]:
        
        return self.eigval, self.eigvec

    def getDecayRates(self, sort : bool = True) -> ndarray:
        """
        TODO: Description
        """
        from numpy import imag, sort

        #if ham is empty
        self.assertInit()
        
        #if eigval is empty
        if not self.decomposed:
            self.eigenDecomposition()
        
        decay_rates = 2 * imag(self.eigval)     #Factor of 2, see Asenjo-Garcia et al.

        if sort:
            decay_rates = sort(decay_rates)

        return decay_rates

    def isDecomposed(self) -> bool:
        return self.decomposed
    
    def isInit(self) -> bool:
        return self.initialized

    """
    SET methods:
    """

    def setN(self, N : int):
        self.N = N

    def setHam(self, ham : ndarray):
        self.hamiltonian = ham

        #flags:
        self.initialized = True
        self.decomposed = False

    """----- Construct standard hamiltonians -----"""

    """
    In this section, the full Hamiltonian is constructed (full 2^N x 2^N space), which is very memory-demanding for larger N. 

    N.B.: Currently incompatible
    """

    def coherence_operators(self, i, j, N):
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

    def H_eff(self, N, w0, D, G):
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
                    H_eff += (-mu_0 * w0**2) * D.trans() * Qobj(G[i,j]) * D * self.coherence_operators(i, j, N)
        return H_eff

    def H(self, N, w0, H_eff):
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

    """
    BLOCK HAMILTONIANS IN BASIS OF {|e_j>} - denoting the single-excitation states of excitation at j'th atom and ground-state elsewhere.
    The Hamiltonian as described in e.g. Asenjo-Garcia commute with the number-operator meaning it conserves excitation-numbers. 
    In other words, the Hamiltonian will be of block-structure nicely divided in 0, 1, 2, 3, ... and so on excitations. 

    Result is arrays of size e.g. (N x N) in the single-excitation case. 
    """

    def block(self, N : int, G : ndarray, n : ndarray):
        """
        Desciption: TODO
        """
        from numpy import zeros
        from scipy.constants import pi

        #Raise flags:
        self.initialized = False
        self.decomposed = False

        self.setN(N)

        #If only one directional vector is passed, e.g. ez, then fill (N x 3) dim array of given direction
        if n.ndim == 1:
            arr = zeros((N,3))
            for i in range(N):
                arr[i] = n
            n = arr
        
        #We want Hamiltonian in basis of {|e_j>}, meaning one excitation on j'th subspace.
        #Do this by computing every matrix element of NxN matrix. It is the block of one excitation in full Hamiltonian
        block = zeros((N,N), dtype=complex)

        #Manually set diagval, see meeting notes 7/3:
        diagval = - 1j/2

        for i in range(N):
            for j in range(N):
                if i == j:
                    block[i, j] = diagval
                else:
                    block[i, j] += -3*pi * n[i].transpose() @ G[i,j] @ n[i]        #Dipoles polarized along z-direction.  

        #Store in internal container.
        self.hamiltonian = block

        #Flag:
        self.initialized = True

    def scalarham(self, N : int, rij : ndarray, w0 : float = 1, dimensionless : bool = True):
        """
        See scalar() under GreensTensor.py
        """
        from GreensTensor import scalar
        from scipy.constants import pi

        #Flags:
        self.initialized = False
        self.decomposed = False

        self.setN(N)

        h = scalar(N, rij, w0, dimensionless)

        #Manually set diagval and multiply with constants, see meeting notes 7/3:
        diagval = - 1j/2
        for i in range(N):
            for j in range(N):
                if i == j:
                    h[i, j] = diagval
                else:
                    h[i, j] = -3 * pi * h[i, j]

        #store in internal
        self.hamiltonian = h
        
        #Lower flag:
        self.initialized = True

    def eigenDecomposition(self):
        """
        TODO: Description
        """
        from numpy.linalg import eig    #import eig, since ham might not be herm

        #Check if ham is empty
        self.assertInit()
        
        eigval, eigvec = eig(self.hamiltonian)

        #Store:
        self.eigval = eigval
        self.eigvec = eigvec

        #Flags:
        self.decomposed = True
        