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
import pylab as plb

plt.style.use('ggplot')
plb.rcParams['font.size'] = 27
plt.rcParams["figure.figsize"] = (16,10)
#%%  PLOTTING CHAIN LENGTH VS GROUND STATE ENERGY
N = 8
bond = 3
tau = 0.01
num_iter = 500
mid_steps = 5
N_sp = [i for i in range(3,13)]
tmp_engs = []
for i in range(0,len(N_sp)):
    print(i)
    res = run_tebd(
        model = 'ising',
        model_params = {'J': 1.0, 'lmda': 0.0},
        N = N_sp[i],
        bond_dim = bond,
        tau = tau,
        num_iter = num_iter,
        mid_steps = mid_steps,
        observables = ["energy"],
        print_to_stdout = False,
        evol_type = "imag",
        st_order = "ST1",
        initial_state = 'random'
        )
    tmp_engs.append(res['energy'][-1])

def analy(N):
    return(-(N**2-N)/N)

plt.figure('Chain')
plt.plot(N_sp, tmp_engs, 'bo', label='TEBD')
plt.plot(N_sp, analy(np.array(N_sp)), label = 'analytic')
plt.xlabel('Chain length N')
plt.ylabel('Ground state energy')
plt.legend()
#plt.savefig('chainlengthvsgsss.pdf', dpi=1000,bbox_inches='tight')
#%%  PLOTTING DIFFERENT TROTTER LEVELS (KEEP SAME INITIAL PARAMETERS, ALSO STATE) (or now trying to do random and average)
N = 8
J = 1.0
lambd = 0.0
bond = 10
tau = 0.5
num_iter = 10
mid_steps = 1
i_temp = 100 #number of iterations to average over

plt.figure('energies')
st1_av = np.zeros(i_temp)
st2_av = np.zeros(i_temp)
st4_av = np.zeros(i_temp)

for i in range(i_temp):
    for st in ['ST1', 'ST2', 'ST4']:
        if st == 'ST1':
            tau =   0.1
            num_iter = i_temp
        if st == 'ST2':
            tau =   0.1
            num_iter = i_temp
        if st == 'ST4':
            tau =   0.1
            num_iter = i_temp
            
        res = run_tebd(
            model = 'ising',
            model_params = {'J': 1.0, 'lmda': 0.0},
            N = N,
            bond_dim = bond,
            tau = tau,
            num_iter = num_iter,
            mid_steps = mid_steps,
            observables = ["energy"],
            print_to_stdout = False,
            evol_type = "imag",
            st_order = st,
            initial_state = 'random'
            )
        if st == 'ST1':
            st1_av += 1/7*res['energy']
        if st == 'ST2':
            st2_av += 1/7*res['energy']
        if st == 'ST4':
            st4_av += 1/7*res['energy']
st1_av = 1/i_temp*st1_av
st2_av = 1/i_temp*st2_av
st4_av = 1/i_temp*st4_av

plt.plot(st1_av, label = 'st1')
plt.plot(st2_av, label = 'st2')
plt.plot(st4_av, label = 'st4')
plt.legend()
#tutorial 3 and 4 of the tensors.net >> canonical form. ?, try run st2 as well and see difference wrt st1
#%% PLOTTING DIFFERENT VALUES FOR TAU (KEEP SAME INITIAL PARAMETERS, ALSO STATE)
N = 8
J = 1.0
lambd = 0.0
bond = 5
tau = 0.01
num_iter = 5*10**4
fig, ax = plt.subplots()
ax.set_xscale('log')
for tau in [0.1, 0.01, 0.001, 0.0001]:
    mid_steps = 5
    if tau == 0.1 or tau == 0.01:
        mid_steps = 1
    res = run_tebd(
        model = 'ising',
        model_params = {'J': 1.0, 'lmda': 0.0},
        N = N,
        bond_dim = bond,
        tau = tau,
        num_iter = num_iter,
        mid_steps = mid_steps,
        observables = ["energy"],
        print_to_stdout = False,
        evol_type = "imag",
        st_order = "ST2",
        initial_state = '11111111'
        )
    ax.plot(res['energy'], label=f'\u03C4 = {tau}')
ax.set_xlim(1, 10**4)
plt.xlabel('Time-step')
plt.ylabel('Energy')
plt.legend()
#%% QUANTUM PHASE TRANSITION
N = 7
J = 1.0
bond = 5
init_state = "0000000"
tau = 0.01
num_iter = 1000
mid_steps = 5
lmd= -np.array([0.0005, 0.01, 0.1,0.2,0.3,0.4, 0.5, 0.75, 1, 1.5, 3, 10, 100, 1000])
M = []
for i in range(0,len(lmd)):
    res = run_tebd(
        model = 'ising',
        model_params = {'J': 1.0, 'lmda': lmd[i]},
        N = N,
        bond_dim = bond,
        tau = tau,
        num_iter = num_iter,
        mid_steps = mid_steps,
        observables = ["magnetization"],
        print_to_stdout = False,
        evol_type = "imag",
        st_order = "ST2",
        initial_state = init_state
        )
    M.append(res['magnetization'][-1])

fig, ax = plt.subplots()

ax.plot(-lmd, M)
ax.plot(-lmd, M, 'bo', label = 'N=7')
plt.xlabel('\u03BB')
plt.ylabel('M')
ax.set_xscale('log')
#plt.savefig('magnetfinalassignment.pdf', dpi=1000,bbox_inches='tight')

#%% THIS IS FOR THE ENTROPY PLOTSS
from matrix_product_states import MatrixProductState

N = 10
J = 1.0
bond = 2
init_state = "ghz"
tau = 0.01
num_iter = 1000
mid_steps = 1
lmd = 1
res = run_tebd(
    model = 'ising',
    model_params = {},
    N = N,
    bond_dim = bond,
    tau = tau,
    num_iter = num_iter,
    mid_steps = mid_steps,
    observables = ["entropy", "energy"],
    print_to_stdout = False,
    evol_type = "real",
    st_order = "ST2",
    initial_state = init_state
    )
surf = plt.pcolormesh(res['entropy'])
plt.grid()
plt.xlabel('site')
plt.ylabel('Time-step')
plt.colorbar(surf, shrink=0.5, aspect=5, label = 'entropy')
#plt.savefig('entropy_unentangledinital.pdf', dpi=1000,bbox_inches='tight')
