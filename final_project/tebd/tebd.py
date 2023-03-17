"""Run TEBD algorithm given various parameters."""

from typing import Optional

import numpy as np
from numpy import linalg as LA
import quimb.tensor as qtn
from scipy.linalg import expm

from tebd.hamiltonian import LocalHamiltonian
from tebd.hamiltonian import Hamiltonian
from tebd.matrix_product_states import MatrixProductState


class TEBD:
    """Class implementing the time-evolving block decimation algorithm."""

    REAL_TIME_EVOLUTION = "real"
    IMAG_TIME_EVOLUTION = "imag"

    ST_ORDER_1 = "ST1"
    ST_ORDER_2 = "ST2"

    def __init__(
            self, mps: MatrixProductState, local_H: LocalHamiltonian, global_H: Hamiltonian, evol_type: str,
            bond_dim: Optional[int] = None, st_order: str = None
    ):
        """Initialize TEBD algorithm with states and Hamiltonian.

        Args:
            mps: Matrix product state object
            H: Hamiltonian object
            evol_type: Type of time evolution (e.g., "real" or "imag")
            bond_dim: Bond dimension to use when updating (if not set, use bond_dim from the underlying MPS)
            st_order: Order of Suzuki-Trotter decomposition (i.e., "ST1" or "ST2")

        """
        if evol_type not in [self.REAL_TIME_EVOLUTION, self.IMAG_TIME_EVOLUTION]:
            raise Exception(f"Evolution type {evol_type} not supported")

        if not st_order:
            st_order = self.ST_ORDER_1
        elif st_order not in [self.ST_ORDER_1, self.ST_ORDER_2]:
            raise Exception(f"Suzuki-Trotter order {st_order} not supported")

        if mps.d != local_H.d:
            raise Exception("Dimension sizes of MPS and Hamiltonian do not agree")

        if mps.N != local_H.N:
            raise Exception("Number of sites of MPS and Hamiltonian do not agree")

        self.mps = mps
        self.local_H = local_H
        self.global_H = global_H
        self.evol_type = evol_type
        self.st_order = st_order

        self.d = mps.d
        self.N = mps.N

        if not bond_dim:
            self.bond_dim = mps.bond_dim
        else:
            self.bond_dim = bond_dim

        # create auxiliary bond matrices
        self.sbonds = [np.ones(self.mps.bond_dim) / np.sqrt(self.mps.bond_dim)] * (self.N - 1)

    def step(self, tau: float):
        """Run one step of TEBD.

        Args:
            tau: Time step

        """
        if self.st_order == self.ST_ORDER_1:
            # apply all gates successively for the same time
            for gate_idx in range(self.N - 1):
                self._apply_gate(gate_idx, tau)
        elif self.st_order == self.ST_ORDER_2:
            # apply odd gates for tau / 2, then even for tau, then odd again for tau / 2
            odd_gate_nums = range(0, self.N - 1, 2)
            even_gate_nums = range(1, self.N - 1, 2)

            for idx in odd_gate_nums:
                self._apply_gate(idx, tau / 2)

            for idx in even_gate_nums:
                self._apply_gate(idx, tau)

            for idx in odd_gate_nums:
                self._apply_gate(idx, tau / 2)

        # renormalize
        self.mps.normalize()

    def compute_energy(self) -> float:
        """Compute energy of the MPS.

        Returns:
            Energy

        """
        rho = self.mps.rho()
        rhoC = rho ^ ...

        inds = tuple([f'k{i}' for i in range(2 * self.N)])

        ham_tensor = qtn.Tensor(self.global_H.hamiltonian, inds=inds, tags=['ham'])
        rho_tensor = qtn.Tensor(rhoC.data, inds=inds, tags=['rho'])

        energy_tensor = ham_tensor & rho_tensor
        energy = energy_tensor ^ ...

        return energy

    def _apply_gate(self, gate_idx: int, tau: float):
        """Apply gate with given index for given time.

        Args:
            gate_idx: Index of gate to apply
            tau: Time step

        """
        two_site_gate = self._gen_gate(self.local_H.hamiltonians[gate_idx], tau)

        if gate_idx == 0:
            # apply left-most gate
            sAB = self._apply_left_gate(
                self.mps.data[gate_idx], self.mps.data[gate_idx + 1], self.sbonds[gate_idx], two_site_gate
            )
        elif gate_idx == self.N - 2:
            # apply right-most gate
            sAB = self._apply_right_gate(
                self.mps.data[gate_idx], self.mps.data[gate_idx + 1], self.sbonds[gate_idx], two_site_gate
            )
        else:
            sAB = self._apply_interior_gate(
                two_site_gate,
                self.mps.data[gate_idx],
                self.mps.data[gate_idx + 1],
                self.sbonds[gate_idx - 1],
                self.sbonds[gate_idx],
                self.sbonds[gate_idx + 1]
            )

        self.sbonds[gate_idx] = sAB

    def _apply_left_gate(
            self, left_site: qtn.Tensor, right_site: qtn.Tensor, central_bond: np.array, gate: np.array
    ) -> np.array:
        """Apply gate to left-most sites.

        Args:
            left_site: Left-most site
            right_site: Site adjacent to left-most site
            central_bond: Bond between two sites
            gate: Gate representing time evolution

        Returns:
            New bond weights

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
        left_site.modify(data=utemp.reshape(self.d, chitemp))
        right_site.modify(data=vhtemp.reshape(chitemp, self.d, right_site.data.shape[2]))
        central_bond = stemp[range(chitemp)] / LA.norm(stemp[range(chitemp)])

        return central_bond

    def _apply_right_gate(
            self, left_site: qtn.Tensor, right_site: qtn.Tensor, central_bond: np.array, gate: np.array
    ) -> np.array:
        """Apply gate to right-most sites.

        Args:
            left_site: Left-most site
            right_site: Site adjacent to left-most site
            central_bond: Bond between two sites
            gate: Gate representing time evolution

        Returns:
            New bond weights

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
        left_site.modify(data=utemp[:, range(chitemp)].reshape(left_site.data.shape[0], self.d, chitemp))
        right_site.modify(data=vhtemp[range(chitemp), :].reshape(chitemp, self.d))
        central_bond = stemp[range(chitemp)] / LA.norm(stemp[range(chitemp)])

        return central_bond

    def _apply_interior_gate(
            self, gate: np.array, left_site: qtn.Tensor, right_site: qtn.Tensor, left_bond: np.array,
            central_bond: np.array, right_bond: np.array, stol=1e-7
    ) -> np.array:
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
            New bond weights

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
        left_site.modify(data=(np.diag(1 / left_bond) @ utemp).reshape(left_site.data.shape[0], self.d, chitemp))
        right_site.modify(data=(vhtemp @ np.diag(1 / right_bond)).reshape(chitemp, self.d, right_site.data.shape[2]))
        central_bond = stemp[range(chitemp)] / LA.norm(stemp[range(chitemp)])

        return central_bond

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


def run_tebd(tebd_obj: TEBD, tau: float, num_iter: int, mid_steps: int):
    """Run the TEBD algorithm for given number of iterations.

    Args:
        tebd_obj: Object representing TEBD algorithm
        tau: Timestep
        num_iter: Number of iterations
        mid_steps: Number of steps between each diagnostic

    """
    energies = []
    wave_functions = []

    for k in range(num_iter):
        if np.mod(k, mid_steps) == 0 or k == num_iter:
            # compute energy
            energy = tebd_obj.compute_energy()
            wave_function = tebd_obj.mps.wave_function()

            print(f"Iteration: {k} of {num_iter}, energy: {energy}")

            energies.append(energy)
            wave_functions.append(wave_function)

        tebd_obj.step(tau)

    wave_functions = np.array(wave_functions)

    return energies, wave_functions
