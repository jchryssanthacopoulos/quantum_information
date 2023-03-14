"""Evolve the matrix product states with local Hamiltonians."""

from numpy import linalg as LA
import numpy as np

import quimb.tensor as qtn


def apply_gate_MPS(gateAB, A, sAB, B, sBA, chi, stol=1e-7):
    """ apply a gate to an MPS across and a A-B link. Truncate the MPS back to
    some desired dimension chi"""

    # ensure singular values are above tolerance threshold
    sBA_trim = sBA * (sBA > stol) + stol * (sBA < stol)

    # contract gate into the MPS, then deompose composite tensor with SVD
    d = A.shape[1]
    chiBA = sBA_trim.shape[0]
    nshape = [d * chiBA, d * chiBA]
    
    sBA_1_tensor = qtn.Tensor(np.diag(sBA), inds=('f0', 'k1'), tags=['sBA', '1'])
    A_tensor = qtn.Tensor(A, inds=('k1', 'k2', 'k3'), tags=['A'])
    sAB_tensor = qtn.Tensor(np.diag(sAB), inds=('k3', 'k4'), tags=['sAB'])
    B_tensor = qtn.Tensor(B, inds=('k4', 'k5', 'k6'), tags=['B'])
    sBA_2_tensor = qtn.Tensor(np.diag(sBA), inds=('k6', 'f3'), tags=['sBA', '2'])
    gate_AB_tensor = qtn.Tensor(gateAB, inds=('f1', 'f2', 'k2', 'k5'), tags=['gateAB'])
    
    TN = sBA_1_tensor & gate_AB_tensor & A_tensor & sAB_tensor & B_tensor & sBA_2_tensor
    TNc = TN ^ ...
    x = TNc.data
    
    utemp, stemp, vhtemp = LA.svd(x.reshape(nshape), full_matrices=False)

    # truncate to reduced dimension
    chitemp = min(chi, len(stemp))
    utemp = utemp[:, range(chitemp)].reshape(sBA_trim.shape[0], d * chitemp)
    vhtemp = vhtemp[range(chitemp), :].reshape(chitemp * d, chiBA)

    # remove environment weights to form new MPS tensors A and B
    A = (np.diag(1 / sBA_trim) @ utemp).reshape(sBA_trim.shape[0], d, chitemp)
    B = (vhtemp @ np.diag(1 / sBA_trim)).reshape(chitemp, d, chiBA)

    # new weights
    sAB = stemp[range(chitemp)] / LA.norm(stemp[range(chitemp)])

    return A, sAB, B
