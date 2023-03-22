"""Class for manipulating matrix product states."""

from typing import List
from typing import Optional

import numpy as np
from numpy import linalg as LA
import quimb.tensor as qtn


class MatrixProductState:
    """Class representing a matrix product state with given number of states."""

    def __init__(self, d: int, N: int, bond_dim: int, states: Optional[List[np.array]] = None):
        """Initialize the MPS.

        Args:
            d: Dimension of each state
            N: Number of states
            bond_dim: Bond dimension between states
            states: Optional states to initialize with

        """
        self.d = d
        self.N = N
        self.bond_dim = bond_dim

        self.data = []
        self.singular_values = []

        if not states:
            states = []
            states.append(np.random.rand(d, bond_dim))
            for i in range(1, N - 1):
                states.append(np.random.rand(bond_dim, d, bond_dim))
            states.append(np.random.rand(bond_dim, d))

        # create left-most state
        self.data.append(qtn.Tensor(states[0], inds=("k0", "i0"), tags=["state 1"]))

        for i in range(1, N - 1):
            self.singular_values.append(np.eye(bond_dim, bond_dim))

            self.data.append(
                qtn.Tensor(self.singular_values[-1], inds=(f"i{2 * (i - 1)}", f"i{2 * i - 1}"), tags=[f"SV {i}"])
            )

            self.data.append(
                qtn.Tensor(states[i], inds=(f"i{2 * i - 1}", f"k{i}", f"i{2 * i}"), tags=[f"state {i + 1}"]))

        # create right-most state
        self.singular_values.append(np.eye(bond_dim, bond_dim))

        self.data.append(
            qtn.Tensor(self.singular_values[-1], inds=(f"i{2 * (N - 2)}", f"i{2 * N - 3}"), tags=[f"SV {N - 1}"])
        )

        self.data.append(qtn.Tensor(states[N - 1], inds=(f"i{2 * N - 3}", f"k{N - 1}"), tags=[f"state {N}"]))

        self.normalize()

    @classmethod
    def init_from_state(cls, state_str: str):
        """Create a specific state like '000' corresponding to provided string.

        Args:
            state_str: State like '00010001' to initialize MPS with

        Returns:
            MPS corresponding to state

        """
        array_map = {
            "0": np.array([1.0, 0.0]),
            "1": np.array([0.0, 1.0])
        }

        N = len(state_str)

        states = []
        for idx, state in enumerate(state_str):
            if idx == 0:
                states.append(array_map[state].reshape(-1, 1))
            elif idx == N - 1:
                states.append(array_map[state].reshape(1, -1))
            else:
                states.append(array_map[state].reshape(1, -1, 1))

        return MatrixProductState(d=2, N=N, bond_dim=1, states=states)

    def wave_function(self):
        """Return wavefunction corresponding to MPS."""
        TNc = qtn.TensorNetwork(self.data) ^ ...
        return TNc.data.reshape((self.d ** self.N))

    def rho(self) -> qtn.TensorNetwork:
        """Compute the density matrix of the MPS."""
        rho_tensors = []

        for idx, ten in enumerate(self.data):
            new_data = ten.data.conj()
            new_inds = self._update_indices(ten.inds)
            rho_tensors.append(qtn.Tensor(new_data, inds=new_inds, tags=[f"state {idx + 1} conj"]))

        for ten in self.data:
            rho_tensors.append(ten.copy())

        return qtn.TensorNetwork(rho_tensors)

    def entropy(self, site: int) -> float:
        """Get entropy between sites to the left of and including given site, and rest of system.

        Args:
            site: End of left bipartition to get entropy for (starts at 1)

        Returns:
            Entropy

        """
        density_matrix = self.rho()

        # identify legs to contract in density matrix
        for i in range(site):
            ten = density_matrix.tensors[2 * i]

            new_inds = ()
            for ind in ten.inds:
                if ind.startswith("k"):
                    new_inds += (f"k{i}",)
                else:
                    new_inds += (ind,)

            density_matrix.tensors[2 * i].modify(inds=new_inds)

        reduced_density_matrix = density_matrix ^ ...

        if isinstance(reduced_density_matrix, float):
            # all indices have been contracted, so density matrix is just a number
            reduced_density_matrix = [[reduced_density_matrix]]
        else:
            # merge indices
            n_reduced = self.d ** (self.N - site)
            reduced_density_matrix = reduced_density_matrix.data.reshape(n_reduced, n_reduced)

        # get eigenvalues
        eigenvalues, _ = LA.eigh(reduced_density_matrix)

        return -sum(w * np.log(w) for w in eigenvalues if w > 0)

    def norm(self):
        """Get the norm of the MPS."""
        tn = qtn.TensorNetwork(self.data)
        norm = tn.H @ tn
        return norm

    def normalize(self):
        """Normalize the MPS."""
        self.multiply(1 / np.sqrt(self.norm()))

    def multiply(self, fac: float):
        """Multiply MPS by given factor.

        Args:
            fac: Factor to multiply by

        """
        self.data[0].modify(data=self.data[0].data * fac)

    def get_state(self, idx: int) -> qtn.Tensor:
        """Get state for site given by index.

        Args:
            idx: Index of state to retrieve

        Returns:
            Tensor corresponding to state

        """
        return self.data[2 * idx]

    def get_sv(self, idx: int) -> qtn.Tensor:
        """Get singular values for given index.

        Args:
            idx: Index of singular values to retrieve

        Returns:
            Tensor corresponding to singular values

        """
        return self.data[2 * idx + 1]

    def _update_indices(self, inds):
        """Update the given indices to correspond to new internal and external legs.

        Args:
            inds (tuple): Indices

        Returns:
            New indices

        """
        new_inds = ()
        for i in inds:
            if i.startswith("k"):
                new_inds += (f"k{int(i[1:]) + self.N}",)
            if i.startswith("i"):
                new_inds += (f"i{int(i[1:]) + 2 * (self.N - 1)}",)
        return new_inds
