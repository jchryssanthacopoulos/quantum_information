"""Run TEBD algorithm given various parameters."""

from typing import Optional

import numpy as np
from scipy.linalg import expm

from tebd.density_matrix import local_density_MPS
from tebd.observables import compute_energy
from tebd.time_evolution import apply_gate_MPS


def run_tebd(
        hamAB: np.ndarray,
        hamBA: np.ndarray,
        A: np.ndarray,
        B: np.ndarray,
        sAB: np.ndarray,
        sBA: np.ndarray,
        chi: int,
        tau: float,
        evotype: Optional[str] = 'imag',
        numiter: Optional[int] = 1000,
        midsteps: Optional[int] = 10,
        E0: Optional[float] = 0.0
):
    """
    Implementation of time evolution (real or imaginary) for MPS with 2-site unit
    cell (A-B), based on TEBD algorithm.
    Args:
    hamAB: nearest neighbor Hamiltonian coupling for A-B sites.
    hamBA: nearest neighbor Hamiltonian coupling for B-A sites.
    A: MPS tensor for A-sites of lattice.
    B: MPS tensor for B-sites of lattice.
    sAB: vector of weights for A-B links.
    sBA: vector of weights for B-A links.
    chi: maximum bond dimension of MPS.
    tau: time-step of evolution.
    evotype: set real (evotype='real') or imaginary (evotype='imag') evolution.
    numiter: number of time-step iterations to take.
    midsteps: number of time-steps between re-orthogonalization of the MPS.
    E0: specify the ground energy (if known).
    Returns:
    np.ndarray: MPS tensor for A-sites;
    np.ndarray: MPS tensor for B-sites;
    np.ndarray: vector sAB of weights for A-B links.
    np.ndarray: vector sBA of weights for B-A links.
    np.ndarray: two-site reduced density matrix rhoAB for A-B sites
    np.ndarray: two-site reduced density matrix rhoAB for B-A sites
    """
    # exponentiate Hamiltonian
    d = A.shape[1]
    if evotype == "real":
        gateAB = expm(1j * tau * hamAB.reshape(d**2, d**2)).reshape(d, d, d, d)
        gateBA = expm(1j * tau * hamBA.reshape(d**2, d**2)).reshape(d, d, d, d)
    elif evotype == "imag":
        gateAB = expm(-tau * hamAB.reshape(d**2, d**2)).reshape(d, d, d, d)
        gateBA = expm(-tau * hamBA.reshape(d**2, d**2)).reshape(d, d, d, d)

    for k in range(numiter + 1):
        if np.mod(k, midsteps) == 0 or (k == numiter):
            """ Compute energy and display """

            # compute 2-site local reduced density matrices
            rhoAB, rhoBA = local_density_MPS(A, sAB, B, sBA)

            energy = compute_energy(hamAB, rhoAB, hamBA, rhoBA)

            chitemp = min(A.shape[0], B.shape[0])
            enDiff = energy - E0

            print('iteration: %d of %d, chi: %d, t-step: %f, energy: %f, '
                'energy error: %e' % (k, numiter, chitemp, tau, energy, enDiff))

        """ Do evolution of MPS through one time-step """
        if k < numiter:
            # apply gate to A-B link
            A, sAB, B = apply_gate_MPS(gateAB, A, sAB, B, sBA, chi)

            # apply gate to B-A link
            B, sBA, A = apply_gate_MPS(gateBA, B, sBA, A, sAB, chi)

    rhoAB, rhoBA = local_density_MPS(A, sAB, B, sBA)
    
    return A, B, sAB, sBA, rhoAB, rhoBA
