"""Use the TEBD algorithm to generate some results and plots."""

from hamiltonian import LocalHamiltonian
from hamiltonian import IsingHamiltonian
from hamiltonian import LocalIsingHamiltonian
from hamiltonian import Hamiltonian
from matrix_product_states import MatrixProductState
from tebd import TEBD
from tebd import run_tebd
import matplotlib.pyplot as plt
import numpy as np
#%%

d = 2
N = 8
J = 1.0
lambd = 0.0
bond = 12
init_state = "11111111"

ham = Hamiltonian(d,N)
ham_isi = IsingHamiltonian(N,J,lambd)

loc_ham = LocalHamiltonian(d,N)
loc_ham_isi = LocalIsingHamiltonian(N,J,lambd)
mps = MatrixProductState(d,N,bond)
for st in ['ST1', 'ST2', 'ST4']:
    mps = mps.init_from_state(init_state)
    tst = TEBD(mps, loc_ham_isi, ham_isi, "imag", bond, st)
    plt.figure('energies')
    res = run_tebd(tst, 0.1, 100, 1)
    plt.plot(res[0], label = st)
    if (st == 'ST2') or (st == 'ST4'):
        diff = np.array(res[0]) - np.array(temp)
        plt.figure(st)
        plt.plot(diff, label = st)
    temp = res[0]
plt.legend()
#tutorial 3 and 4 of the tensors.net >> canonical form. ?, try run st2 as well and see difference wrt st1
#%%
fig, ax = plt.subplots()
ax.set_xscale('log')
mps = MatrixProductState(d,N,bond)

for tau in [0.1, 0.01, 0.001, 0.0001]:
    mps = mps.init_from_state(init_state)
    tst = TEBD(mps, loc_ham_isi, ham_isi, "imag", bond, 'ST2')
    res = run_tebd(tst, tau, 13**4, 1)
    ax.plot(res[0], label = ('$\tau$ =', tau))
    
plt.legend()
#%%
d = 2
N = 7
J = 1.0
bond = 5
init_state = "0000000"

ham = Hamiltonian(d,N)

loc_ham = LocalHamiltonian(d,N)
mps = MatrixProductState(d,N,bond)

lmd= -np.array([0.0005, 0.01, 0.1, 0.5, 0.75, 1, 1.5, 3, 10, 100, 1000])
M = []
for i in range(0,len(lmd)):
    ham_isi = IsingHamiltonian(N,J,lmd[i])
    loc_ham_isi = LocalIsingHamiltonian(N,J,lmd[i])

    mps = mps.init_from_state(init_state)
    tst = TEBD(mps, loc_ham_isi, ham_isi, "imag", bond, 'ST1')
    res = run_tebd(tst, 0.1, 100, 1)
    M.append(res[1][-1])

fig, ax = plt.subplots()

ax.plot(-lmd, M)
ax.plot(-lmd, M, 'bo', label = 'N=7')
plt.xlabel('\u03BB')
plt.ylabel('M')
ax.set_xscale('log')
#plt.savefig('magnetfinalassignment.pdf', dpi=1000,bbox_inches='tight')
