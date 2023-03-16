"""Class describing Hamiltonians of various types of systems."""

import numpy as np

import quimb as qu


class Hamiltonian:
    """Hamiltonian base class."""

    def __init__(self, d: int, N: int):
        self.d = d
        self.N = N
        self.hamiltonians = np.zeros((N - 1, d ** 2, d ** 2))


class IsingHamiltonian(Hamiltonian):
    """Create local Hamiltonians for Ising chain."""

    def __init__(self, N: int, lmda: float):
        """Create Hamiltonians for given number of sites.

        Args:
            N: Number of sites
            lmda: Coupling to external magnetic field

        """
        super(IsingHamiltonian, self).__init__(2, N)

        hamiltonian_two_site = (
            np.kron(qu.pauli("X"), qu.pauli("X")) +
            lmda * np.kron(qu.pauli("Z"), np.eye(2)) +
            lmda * np.kron(np.eye(2), qu.pauli("Z"))
        )

        self.hamiltonians = np.repeat(hamiltonian_two_site[None, :], N - 1, axis=0)
