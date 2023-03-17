"""Class describing Hamiltonians of various types of systems."""

import numpy as np
import quimb as qu


class LocalHamiltonian:
    """Base class for representing local Hamiltonians."""

    def __init__(self, d: int, N: int):
        """Initialize local Hamiltonian.

        Args:
            d: Number of dimensions
            N: Number of sites

        """
        self.d = d
        self.N = N
        self.hamiltonians = np.zeros((N - 1, d ** 2, d ** 2))


class LocalIsingHamiltonian(LocalHamiltonian):
    """Create local Hamiltonians for Ising chain."""

    def __init__(self, N: int, lmda: float):
        """Create Hamiltonians for given number of sites.

        Args:
            N: Number of sites
            lmda: Coupling to external magnetic field

        """
        super(LocalIsingHamiltonian, self).__init__(2, N)

        hamiltonian_two_site = (
            np.kron(qu.pauli("X"), qu.pauli("X")) +
            lmda * np.kron(qu.pauli("Z"), np.eye(2)) +
            lmda * np.kron(np.eye(2), qu.pauli("Z"))
        )

        self.hamiltonians = np.repeat(hamiltonian_two_site[None, :], N - 1, axis=0)


class Hamiltonian:
    """Base class for representing general Hamiltonians."""

    def __init__(self, d: int, N: int):
        """Initialize Hamiltonian.

        Args:
            d: Number of dimensions
            N: Number of sites

        """
        self.d = d
        self.N = N
        self.hamiltonian = np.zeros((d ** N, d ** N))


class IsingHamiltonian(Hamiltonian):
    """Create Hamiltonian for Ising chain."""

    def __init__(self, N: int, lmda: float):
        """Create Hamiltonian for given number of sites.

        Args:
            N: Number of sites
            lmda: Coupling to external magnetic field

        """
        super(IsingHamiltonian, self).__init__(2, N)

        H_0 = self.non_interacting_hamiltonian()
        H_int = self.interacting_hamiltonian()

        self.hamiltonian = lmda * H_0 + H_int

    def non_interacting_hamiltonian(self) -> np.array:
        """Get non-interacting part of the Hamiltonian.

        Returns:
            Multi-dimensional array representing non-interacting Hamiltonian

        """
        dim = 2 ** self.N

        H_0 = np.zeros((dim, dim), dtype='complex128')

        for i in range(self.N):
            I1 = self.identity(i)
            I2 = self.identity(self.N - i - 1)
            prod1 = np.kron(I1, qu.pauli("Z"))
            H_0 += np.kron(prod1, I2)

        H_0 = H_0.reshape((self.d,) * 2 * self.N)

        return H_0

    def interacting_hamiltonian(self) -> np.array:
        """Get interacting part of the Hamiltonian.

        Returns:
            Multi-dimensional array representing interacting Hamiltonian

        """
        dim = 2 ** self.N

        H_int = np.zeros((dim, dim), dtype='complex128')

        for i in range(self.N - 1):
            I1 = self.identity(i)
            I2 = self.identity(self.N - i - 2)
            prod1 = np.kron(I1, qu.pauli("X"))
            prod2 = np.kron(prod1, qu.pauli("X"))

            H_int += np.kron(prod2, I2)

        H_int = H_int.reshape((self.d,) * 2 * self.N)

        return H_int

    def identity(self, n: int) -> np.array:
        """Get identity matrix for n-qubit state.

        Args:
            n: Number of qubits

        Returns:
            Identity matrix for n qubits

        """
        return np.eye(2 ** n)
