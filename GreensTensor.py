"""
Author: tzs820, s240670, Emil Henningsen.

GREEN'S TENSOR (in free space)
In this section, the Green's tensor for the i'th dipole w.r.t. the j'th dipole is calculated.
The assumptions include:
 * Markovian regime: I.e. the bandwidth of atomic transition, w0, is small and the retardation of the field is negligible.
 * Free space: I.e. the Green's tensor is calculated in free space, where the dipoles are tightly trapped, i.e. d < lambda_0.
 * The stroung coupling regime is avoided: I.e. the electromagnetic field is not strongly coupled to the atomic transition (CQED).
 * The dipoles are point-like particles: I.e. r_i and r_j can be treated as fixed position vectors. 
"""

def G_0(r, w):
    """
    Green's tensor in free space. As given by eq. (6) of Asenjo-Garcia et al.

    ---Parameters:
    r: vector (3) of floats, displacement between the dipoles.
    w: float, frequency of the atomic transition.

    ---Return:
    G: array (3 x 3) of complex floats, Green's tensor in free space between two dipoles.
    """
    from numpy.linalg import norm
    from numpy import e, outer, pi
    from scipy.constants import c
    from qutip import qeye
    
    k0 = w/c
    r_norm = norm(r)
    G = ((e**(k0*r_norm*complex(0,1)))/(4 * pi * k0**2 * r_norm**3)) * (
        (k0**2 * r_norm**2 + k0 * r_norm * complex(0, 1) - 1)*qeye(3) 
        + (-k0**2 * r_norm**2 - 3 * k0 * r_norm * complex(0, 1) + 3)*outer(r, r)/(r_norm**2)
        )
    return G

def fill_G(N, rij, w0, type={"free_space"}):
    """
    Fill array of chosen Green's Tensor's for a given lattice with displacement vectors, rij. 

    ---Parameters:
    N: Integer number of dipoles.
    rij: array (N x N x 3) of floats, displacement vectors between dipoles in lattice.
    w0: float, dipole transition frequency.
    type: string, choose which Green's Tensor to compute.

    ---Return.
    G: array (N x N x 3 x 3) of complex floats, the Green's Tensor in free space for given lattice. 
    """
    from numpy import zeros, eye

    G = zeros((N, N, 3, 3), dtype=complex)                        #Array of Green's tensors
    for i in range(N):                                              #Fill the array of Green's tensors
        for j in range(i, N):
            if i == j:
                G[i, j] = eye(3, dtype=complex)              #Set diagonal values to 1 to avoid divergence.
            else:
                if type == "free_space":                          #Choose Green's Tensor by given optional argument.
                    G[i, j] = G_0(rij[i, j], w0)
                G[j, i] = G[i, j]                                 #By construction, Green's tensor is symmetric.
    return G
