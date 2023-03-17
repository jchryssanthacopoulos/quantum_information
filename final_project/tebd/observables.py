"""Compute various observables."""

import numpy as np
import quimb.tensor as qtn


def compute_energy(hamAB: np.array, rhoAB: np.array, hamBA: np.array, rhoBA: np.array) -> float:
    """Calculate the energy for the given Hamiltonian and density matrix.

    Args:
        hamAB: Hamiltonian for A-B link
        rhoAB: Density matrix for A-B link
        hamBA: Hamiltonian for B-A link
        rhoBA: Density matrix for B-A link

    Returns:
        Average energy

    """
    hamAB_tensor = qtn.Tensor(hamAB, inds=('k1', 'k2', 'k3', 'k4'), tags=['hamAB'])
    rhoAB_tensor = qtn.Tensor(rhoAB, inds=('k1', 'k2', 'k3', 'k4'), tags=['rhoAB'])
    energyAB_tensor = hamAB_tensor & rhoAB_tensor
    energyAB = energyAB_tensor ^ ...

    hamBA_tensor = qtn.Tensor(hamBA, inds=('k1', 'k2', 'k3', 'k4'), tags=['hamBA'])
    rhoBA_tensor = qtn.Tensor(rhoBA, inds=('k1', 'k2', 'k3', 'k4'), tags=['rhoBA'])
    energyBA_tensor = hamBA_tensor & rhoBA_tensor
    energyBA = energyBA_tensor ^ ...

    energy = 0.5 * (energyAB + energyBA)

    return energy
