"""Class for manipulating matrix product states."""

import numpy as np

import quimb.tensor as qtn


class MatrixProductState:
    """Class representing a matrix product state with given number of states."""

    def __init__(self, d: int, N: int, bond_dim: int):
        """Initialize the MPS.

        Args:
            d: Dimension of each state
            N: Number of states
            bond_dim: Bond dimension between states

        """
        if N < 3:
            raise Exception("Number of states must be at least 3")

        self.d = d
        self.N = N
        self.bond_dim = bond_dim

        self.data = []

        # create left-most tensor
        state = np.random.rand(d, bond_dim)
        self.data.append(qtn.Tensor(state, inds=("k0", "i0"), tags=["state 1"]))

        for i in range(1, N - 1):
            state = np.random.rand(bond_dim, d, bond_dim)
            self.data.append(qtn.Tensor(state, inds=(f"i{i - 1}", f"k{i}", f"i{i}"), tags=[f"state {i + 1}"]))

        # create right-most state
        state = np.random.rand(bond_dim, d)
        self.data.append(qtn.Tensor(state, inds=(f"i{N - 2}", f"k{N - 1}"), tags=[f"state {N}"]))

    def rho(self):
        rho_tensors = []

        for idx, ten in enumerate(self.data):
            new_data = ten.data.conj()
            new_inds = self._update_indices(ten.inds)
            rho_tensors.append(qtn.Tensor(new_data, inds=new_inds, tags=[f"state {idx + 1} conj"]))

        for ten in self.data:
            rho_tensors.append(ten.copy())

        return qtn.TensorNetwork(rho_tensors)

    def _update_indices(self, inds):
        new_inds = ()
        for i in inds:
            if i.startswith("k"):
                new_inds += (f"k{int(i[1:]) + self.N}",)
            if i.startswith("i"):
                new_inds += (f"i{int(i[1:]) + self.N - 1}",)
        return new_inds
