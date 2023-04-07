"""Run TEBD algorithm given various parameters."""

from typing import Dict
from typing import List
from typing import Optional

import numpy as np
from numpy import linalg as LA
import quimb.tensor as qtn
from scipy.linalg import expm

from tebd.hamiltonian import LocalHamiltonian
from tebd.hamiltonian import LocalHeisenbergHamiltonian
from tebd.hamiltonian import LocalIsingHamiltonian
from tebd.hamiltonian import Hamiltonian
from tebd.hamiltonian import HeisenbergHamiltonian
from tebd.hamiltonian import IsingHamiltonian
from tebd.matrix_product_states import MatrixProductState


class TEBD:
    """Class implementing the time-evolving block decimation algorithm."""

    REAL_TIME_EVOLUTION = "real"
    IMAG_TIME_EVOLUTION = "imag"

    ST_ORDER_1 = "ST1"
    ST_ORDER_2 = "ST2"
    ST_ORDER_4 = "ST4"

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
            st_order: Order of Suzuki-Trotter decomposition (i.e., "ST1", "ST2", "ST4")

        """
        if evol_type not in [self.REAL_TIME_EVOLUTION, self.IMAG_TIME_EVOLUTION]:
            raise Exception(f"Evolution type {evol_type} not supported")

        if not st_order:
            st_order = self.ST_ORDER_1
        elif st_order not in [self.ST_ORDER_1, self.ST_ORDER_2, self.ST_ORDER_4]:
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
        elif self.st_order == self.ST_ORDER_4:
            tau_1 = (1 / (4 - 4 ** (1 / 3))) * tau
            tau_2 = (1 - 4 * tau_1) * tau
            odd_gate_nums = range(0, self.N - 1, 2)
            even_gate_nums = range(1, self.N - 1, 2)

            for i in range(0, 2):
                for idx in odd_gate_nums:
                    self._apply_gate(idx, tau_1 / 2)

                for idx in even_gate_nums:
                    self._apply_gate(idx, tau_1)

                for idx in odd_gate_nums:
                    self._apply_gate(idx, tau_1 / 2)

            for idx in odd_gate_nums:
                self._apply_gate(idx, tau_2 / 2)

            for idx in even_gate_nums:
                self._apply_gate(idx, tau_2)

            for idx in odd_gate_nums:
                self._apply_gate(idx, tau_2 / 2)

            for i in range(0, 2):
                for idx in odd_gate_nums:
                    self._apply_gate(idx, tau_1 / 2)

                for idx in even_gate_nums:
                    self._apply_gate(idx, tau_1)

                for idx in odd_gate_nums:
                    self._apply_gate(idx, tau_1 / 2)

        # renormalize
        self.mps.normalize()

    def sweep(
            self, tau: float, num_iter: int, mid_steps: int, observables: Optional[List[str]] = None,
            print_to_stdout: Optional[bool] = True
    ) -> Dict:
        """Run the TEBD algorithm for given number of iterations.

        Args:
            tau: Timestep
            num_iter: Number of iterations
            mid_steps: Number of steps between each diagnostic
            observables: List of observables to compute at each midstep
            print_to_stdout: Whether to print diagnostic information to screen

        Returns:
            Dictionary of observables at each midstep

        """
        observables_at_midsteps = {k: [] for k in observables or []}

        directions = ["x", "y", "z"]

        for k in range(num_iter):
            if np.mod(k, mid_steps) == 0:
                # compute observables
                if "energy" in observables_at_midsteps:
                    energy = self.compute_energy()
                    observables_at_midsteps["energy"].append(energy)

                if "entropy" in observables_at_midsteps:
                    entropy = [self.mps.entropy(s) for s in range(1, self.N)]
                    observables_at_midsteps["entropy"].append(entropy)

                if "wavefunction" in observables_at_midsteps:
                    observables_at_midsteps["wavefunction"].append(self.mps.wave_function())

                for direction in directions:
                    if f"{direction} magnetization" in observables_at_midsteps:
                        magnetization = [self.mps.magnetization(n, direction) for n in range(1, self.N + 1)]
                        observables_at_midsteps[f"{direction} magnetization"].append(magnetization)

                if print_to_stdout and "energy" in observables_at_midsteps:
                    print(f"Iteration: {k} of {num_iter}, energy: {energy}")

            self.step(tau)

        observables_at_midsteps = {k: np.array(v) for k, v in observables_at_midsteps.items()}

        return observables_at_midsteps

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
            self._apply_left_gate(
                self.mps.get_state(gate_idx), self.mps.get_state(gate_idx + 1), self.mps.get_sv(gate_idx), two_site_gate
            )
            return

        if gate_idx == self.N - 2:
            # apply right-most gate
            self._apply_right_gate(
                self.mps.get_state(gate_idx), self.mps.get_state(gate_idx + 1), self.mps.get_sv(gate_idx), two_site_gate
            )
            return

        self._apply_interior_gate(
            two_site_gate,
            self.mps.get_state(gate_idx),
            self.mps.get_state(gate_idx + 1),
            self.mps.get_sv(gate_idx - 1),
            self.mps.get_sv(gate_idx),
            self.mps.get_sv(gate_idx + 1)
        )

    def _apply_left_gate(self, left_site: qtn.Tensor, right_site: qtn.Tensor, central_bond: qtn.Tensor, gate: np.array):
        """Apply gate to left-most sites.

        Args:
            left_site: Left-most site
            right_site: Site adjacent to left-most site
            central_bond: Bond between two sites
            gate: Gate representing time evolution

        """
        left_site_T = qtn.Tensor(left_site.data, inds=('k1', 'k2'), tags=['left site'])
        central_bond_T = qtn.Tensor(central_bond.data, inds=('k2', 'k3'), tags=['central bond'])
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
        central_bond.modify(data=np.diag(stemp[range(chitemp)] / LA.norm(stemp[range(chitemp)])))

    def _apply_right_gate(
            self, left_site: qtn.Tensor, right_site: qtn.Tensor, central_bond: qtn.Tensor, gate: np.array
    ):
        """Apply gate to right-most sites.

        Args:
            left_site: Left-most site
            right_site: Site adjacent to left-most site
            central_bond: Bond between two sites
            gate: Gate representing time evolution

        """
        left_site_T = qtn.Tensor(left_site.data, inds=('k1', 'k2', 'k3'), tags=['left site'])
        central_bond_T = qtn.Tensor(central_bond.data, inds=('k3', 'k4'), tags=['central bond'])
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
        central_bond.modify(data=np.diag(stemp[range(chitemp)] / LA.norm(stemp[range(chitemp)])))

    def _apply_interior_gate(
            self, gate: np.array, left_site: qtn.Tensor, right_site: qtn.Tensor, left_bond: qtn.Tensor,
            central_bond: qtn.Tensor, right_bond: np.array, stol=1e-7
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

        """
        # ensure singular values are above tolerance threshold
        left_bond_data = np.diagonal(left_bond.data)
        left_bond_data = np.diag(left_bond_data * (left_bond_data > stol) + stol * (left_bond_data < stol))

        right_bond_data = np.diagonal(right_bond.data)
        right_bond_data = np.diag(right_bond_data * (right_bond_data > stol) + stol * (right_bond_data < stol))

        left_bond_T = qtn.Tensor(left_bond_data, inds=('f0', 'k1'), tags=['left bond'])
        left_site_T = qtn.Tensor(left_site.data, inds=('k1', 'k2', 'k3'), tags=['left site'])
        central_bond_T = qtn.Tensor(central_bond.data, inds=('k3', 'k4'), tags=['central bond'])
        right_site_T = qtn.Tensor(right_site.data, inds=('k4', 'k5', 'k6'), tags=['right site'])
        right_bond_T = qtn.Tensor(right_bond_data, inds=('k6', 'f3'), tags=['right bond'])
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

        # remove environment weights to form new MPS tensors
        left_site.modify(data=(LA.inv(left_bond_data) @ utemp).reshape(left_site.data.shape[0], self.d, chitemp))
        right_site.modify(data=(vhtemp @ LA.inv(right_bond_data)).reshape(chitemp, self.d, right_site.data.shape[2]))
        central_bond.modify(data=np.diag(stemp[range(chitemp)] / LA.norm(stemp[range(chitemp)])))

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
        model: str,
        model_params: Dict,
        N: int,
        bond_dim: Optional[int] = 2,
        tau: Optional[float] = 0.01,
        num_iter: Optional[int] = 500,
        mid_steps: Optional[int] = 10,
        observables: Optional[List[str]] = None,
        print_to_stdout: Optional[bool] = True,
        evol_type: Optional[str] = "imag",
        st_order: Optional[str] = "ST1",
        initial_state: Optional[str] = "random",
        rng_seed: Optional[int] = 0
) -> Dict:
    """Run TEBD for provided model.

    Args:
        model: Model to run
        model_params: Set of model parameters
        N: Number of sites
        bond_dim: Bond dimension
        tau: Timestep
        num_iter: Number of iterations
        mid_steps: Number of steps between each diagnostic
        observables: List of observables to compute at each midstep
        print_to_stdout: Whether to print diagnostic information to screen
        evol_type: Type of time evolution (e.g., "real" or "imag")
        st_order: Order of Suzuki-Trotter decomposition (i.e., "ST1" or "ST2")
        initial_state: What to set as the initial state (either "random" or a string of 1s and 0s of length N)
        rng_seed: Number to seed random number generator

    Returns:
        Dictionary of observables at each midstep

    """
    d = 2

    if initial_state == "random":
        MPS = MatrixProductState(d=d, N=N, bond_dim=bond_dim, rng_seed=rng_seed)
    else:
        MPS = MatrixProductState.init_from_state(initial_state)

    # create Hamiltonians
    if model == "ising":
        loc_ham = LocalIsingHamiltonian(N=N, **model_params)
        glob_ham = IsingHamiltonian(N=N, **model_params)
    elif model == "heisenberg":
        loc_ham = LocalHeisenbergHamiltonian(N=N, **model_params)
        glob_ham = HeisenbergHamiltonian(N=N, **model_params)
    else:
        raise Exception(f"Model {model} not supported")

    # create TEBD object
    tebd_obj = TEBD(MPS, loc_ham, glob_ham, bond_dim=bond_dim, evol_type=evol_type, st_order=st_order)

    # run algorithm
    observables_at_midsteps = tebd_obj.sweep(tau, num_iter, mid_steps, observables, print_to_stdout)

    return observables_at_midsteps
