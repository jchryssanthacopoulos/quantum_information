"""Compute the density matrices associated with matrix product states."""

import numpy as np

import quimb.tensor as qtn


def local_density_MPS(A, sAB, B, sBA):
    """Compute the local reduced density matrix for two sites A and B."""

    # recast singular weights into a matrix
    mAB = np.diag(sAB)
    mBA = np.diag(sBA)

    # contract MPS for local reduced density matrix (A-B)
    sBA_1_tensor = qtn.Tensor(np.diag(sBA**2), inds=('k3', 'k4'), tags=['sBA', '1'])
    A_tensor = qtn.Tensor(A, inds=('k3', 'f3', 'k1'), tags=['A'])
    A_conj_tensor = qtn.Tensor(A.conj(), inds=('k4', 'f1', 'k2'), tags=['A_conj'])
    mAB_1_tensor = qtn.Tensor(mAB, inds=('k1', 'k7'), tags=['mAB', '1'])
    mAB_2_tensor = qtn.Tensor(mAB, inds=('k2', 'k8'), tags=['mAB', '2'])
    B_tensor = qtn.Tensor(B, inds=('k7', 'f4', 'k5'), tags=['B'])
    B_conj_tensor = qtn.Tensor(B.conj(), inds=('k8', 'f2', 'k6'), tags=['B_conj'])
    sBA_2_tensor = qtn.Tensor(np.diag(sBA**2), inds=('k5', 'k6'), tags=['sBA', '2'])
    
    TN_rhoAB = (
        sBA_1_tensor & A_tensor & A_conj_tensor & mAB_1_tensor & mAB_2_tensor & B_tensor & B_conj_tensor & sBA_2_tensor
    )
    TN_rhoAB = TN_rhoAB ^ ...
    rhoAB = TN_rhoAB.transpose('f1', 'f2', 'f3', 'f4').data

    sAB_1_tensor = qtn.Tensor(np.diag(sAB**2), inds=('k3', 'k4'), tags=['sAB', '1'])
    B_tensor = qtn.Tensor(B, inds=('k3', 'f3', 'k1'), tags=['B'])
    B_conj_tensor = qtn.Tensor(B.conj(), inds=('k4', 'f1', 'k2'), tags=['B_conj'])
    mBA_1_tensor = qtn.Tensor(mBA, inds=('k1', 'k7'), tags=['mBA', '1'])
    mBA_2_tensor = qtn.Tensor(mBA, inds=('k2', 'k8'), tags=['mBA', '2'])
    A_tensor = qtn.Tensor(A, inds=('k7', 'f4', 'k5'), tags=['A'])
    A_conj_tensor = qtn.Tensor(A.conj(), inds=('k8', 'f2', 'k6'), tags=['A_conj'])
    sAB_2_tensor = qtn.Tensor(np.diag(sAB**2), inds=('k5', 'k6'), tags=['sAB', '2'])
    
    TN_rhoBA = (
        sAB_1_tensor & B_tensor & B_conj_tensor & mBA_1_tensor & mBA_2_tensor & A_tensor & A_conj_tensor & sAB_2_tensor
    )
    TN_rhoBA = TN_rhoBA ^ ...
    rhoBA = TN_rhoBA.transpose('f1', 'f2', 'f3', 'f4').data

    return rhoAB, rhoBA
