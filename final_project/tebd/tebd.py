"""Run TEBD algorithm given various parameters."""

from typing import Optional

import numpy as np
from numpy import linalg as LA
import quimb.tensor as qtn
from scipy.linalg import expm

from tebd.density_matrix import local_density_MPS
from tebd.observables import compute_energy
from tebd.time_evolution import apply_gate_MPS

from tebd.matrix_product_states import MatrixProductState
from tebd.hamiltonian import Hamiltonian


class TEBD:
    """Class implementing the time-evolving block decimation algorithm."""

    REAL_TIME_EVOLUTION = "real"
    IMAG_TIME_EVOLUTION = "imag"

    def __init__(self, mps: MatrixProductState, H: Hamiltonian, evol_type: str):
        """Initialize TEBD algorithm with states and Hamiltonian.

        Args:
            mps: Matrix product state object
            H: Hamiltonian object
            evol_type: Type of time evolution (e.g., "real" or "imag")

        """
        if evol_type not in [self.REAL_TIME_EVOLUTION, self.IMAG_TIME_EVOLUTION]:
            raise Exception(f"Evolution type {evol_type} not supported")

        if mps.d != H.d:
            raise Exception("Dimension sizes of MPS and Hamiltonian do not agree")

        if mps.N != H.N:
            raise Exception("Number of sites of MPS and Hamiltonian do not agree")

        self.mps = mps
        self.H = H
        self.evol_type = evol_type

        self.d = mps.d
        self.N = mps.N
        self.bond_dim = mps.bond_dim

        # create auxiliary bond matrices
        self.sbonds = [np.ones(self.mps.bond_dim) / np.sqrt(self.mps.bond_dim)] * (self.N - 1)

    def step(self, tau: float):
        """Run one step of TEBD.

        Args:
            tau: Time step

        """
        # apply left-most gate
        two_site_gate = self._gen_gate(self.H.hamiltonians[0], tau)
        A, B, sAB = self._apply_left_gate(self.mps.data[0], self.mps.data[1], self.sbonds[0], two_site_gate)

        self.mps.data[0].modify(data=A)
        self.mps.data[1].modify(data=B)
        self.sbonds[0] = sAB

        for i in range(1, self.N - 2):
            two_site_gate = self._gen_gate(self.H.hamiltonians[i], tau)

            A, B, sAB = self._apply_gate(
                two_site_gate,
                self.mps.data[i],
                self.mps.data[i + 1],
                self.sbonds[i - 1],
                self.sbonds[i],
                self.sbonds[i + 1]
            )

            self.mps.data[i].modify(data=A)
            self.mps.data[i + 1].modify(data=B)
            self.sbonds[i] = sAB
        
        # apply right-most gate
        two_site_gate = self._gen_gate(self.H.hamiltonians[self.N - 2], tau)
        A, B, sAB = self._apply_right_gate(
            self.mps.data[self.N - 2], self.mps.data[self.N - 1], self.sbonds[self.N - 2], two_site_gate
        )

        self.mps.data[self.N - 2].modify(data=A)
        self.mps.data[self.N - 1].modify(data=B)
        self.sbonds[self.N - 2] = sAB

    def _apply_left_gate(self, left_site, right_site, central_bond, gate):
        """Apply gate to left-most sites.
        
        Args:
            left_site: Left-most site
            right_site: Site adjacent to left-most site
            central_bond: Bond between two sites
            gate: Gate representing time evolution

        Returns:
            Tuple of new left, right, and bond

        """
        left_site_T = qtn.Tensor(left_site.data, inds=('k1', 'k2'), tags=['left site'])
        central_bond_T = qtn.Tensor(np.diag(central_bond), inds=('k2', 'k3'), tags=['central bond'])
        right_site_T = qtn.Tensor(right_site.data, inds=('k3', 'k4', 'k5'), tags=['right site'])
        gate_T = qtn.Tensor(gate, inds=('f0', 'f1', 'k1', 'k4'), tags=['gate'])

        TN = gate_T & left_site_T & central_bond_T & right_site_T
        TNc = TN ^ ...

        nshape = [self.d, self.d * right_site.data.shape[2]]
        utemp, stemp, vhtemp = LA.svd(TNc.data.reshape(nshape), full_matrices=False)

        # truncate to reduced dimension
        chitemp = min(self.bond_dim, len(stemp))
        left_site = utemp.reshape(self.d, chitemp)
        right_site = vhtemp.reshape(chitemp, self.d, right_site.data.shape[2])

        # new weights
        central_bond = stemp[range(chitemp)] / LA.norm(stemp[range(chitemp)])

        return left_site, right_site, central_bond

    def _apply_right_gate(self, left_site, right_site, central_bond, gate):
        """Apply gate to right-most sites.
        
        Args:
            left_site: Left-most site
            right_site: Site adjacent to left-most site
            central_bond: Bond between two sites
            gate: Gate representing time evolution

        Returns:
            Tuple of new left, right, and bond

        """
        left_site_T = qtn.Tensor(left_site.data, inds=('k1', 'k2', 'k3'), tags=['left site'])
        central_bond_T = qtn.Tensor(np.diag(central_bond), inds=('k3', 'k4'), tags=['central bond'])
        right_site_T = qtn.Tensor(right_site.data, inds=('k4', 'k5'), tags=['right site'])
        gate_T = qtn.Tensor(gate, inds=('f0', 'f1', 'k2', 'k5'), tags=['gate'])

        TN = left_site_T & gate_T & central_bond_T & right_site_T
        TNc = TN ^ ...

        nshape = [left_site.data.shape[0] * self.d, self.d]
        utemp, stemp, vhtemp = LA.svd(TNc.data.reshape(nshape), full_matrices=False)

        # truncate to reduced dimension
        chitemp = min(self.bond_dim, len(stemp))
        left_site = utemp[:, range(chitemp)].reshape(left_site.data.shape[0], self.d, chitemp)
        right_site = vhtemp[range(chitemp), :].reshape(chitemp, self.d)

        # new weights
        central_bond = stemp[range(chitemp)] / LA.norm(stemp[range(chitemp)])

        return left_site, right_site, central_bond

    def _apply_gate(
            self, gate: np.array, left_site: np.array, right_site: np.array, left_bond: np.array,
            central_bond: np.array, right_bond: np.array, stol=1e-7
    ):
        """Apply gate to two interior sites.

        Args:
            left_site: Left-most site
            right_site: Site adjacent to left-most site
            left_bond: Left to the left of left site
            central_bond: Bond between two sites
            right_bond: Bond to the right of right site
            gate: Gate representing time evolution
            stol: Threshold for singular values

        Returns:
            Tuple of new left, right, and central bond

        """
        # ensure singular values are above tolerance threshold
        left_bond = left_bond * (left_bond > stol) + stol * (left_bond < stol)
        right_bond = right_bond * (right_bond > stol) + stol * (right_bond < stol)

        left_bond_T = qtn.Tensor(np.diag(left_bond), inds=('f0', 'k1'), tags=['left bond'])
        left_site_T = qtn.Tensor(left_site.data, inds=('k1', 'k2', 'k3'), tags=['left site'])
        central_bond_T = qtn.Tensor(np.diag(central_bond), inds=('k3', 'k4'), tags=['central bond'])
        right_site_T = qtn.Tensor(right_site.data, inds=('k4', 'k5', 'k6'), tags=['right site'])
        right_bond_T = qtn.Tensor(np.diag(right_bond), inds=('k6', 'f3'), tags=['right bond'])
        gate_T = qtn.Tensor(gate, inds=('f1', 'f2', 'k2', 'k5'), tags=['gate'])

        # contract with gate
        TN = left_bond_T & gate_T & left_site_T & central_bond_T & right_site_T & right_bond_T
        TNc = TN ^ ...

        # perform SVD
        nshape = [self.d * left_site.data.shape[0], self.d * right_site.data.shape[2]]
        utemp, stemp, vhtemp = LA.svd(TNc.data.reshape(nshape), full_matrices=False)

        # truncate to reduced dimension
        chitemp = min(self.bond_dim, len(stemp))
        utemp = utemp[:, range(chitemp)].reshape(left_site.data.shape[0], self.d * chitemp)
        vhtemp = vhtemp[range(chitemp), :].reshape(chitemp * self.d, right_site.data.shape[2])

        # remove environment weights to form new MPS tensors A and B
        left_site = (np.diag(1 / left_bond) @ utemp).reshape(left_site.data.shape[0], self.d, chitemp)
        right_site = (vhtemp @ np.diag(1 / right_bond)).reshape(chitemp, self.d, right_site.data.shape[2])

        # new weights
        central_bond = stemp[range(chitemp)] / LA.norm(stemp[range(chitemp)])

        return left_site, right_site, central_bond

    def _gen_gate(self, hamiltonian: np.array, tau: float) -> np.array:
        """Generate gate for given Hamiltonian and time step.

        Args:
            hamiltonian: Matrix representing Hamiltonian
            tau: Time step

        Returns:
            Matrix representing time evolution operator

        """
        if self.evol_type == self.REAL_TIME_EVOLUTION:
            return expm(1j * tau * hamiltonian).reshape(self.d, self.d, self.d, self.d)
        return expm(-tau * hamiltonian).reshape(self.d, self.d, self.d, self.d)


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
