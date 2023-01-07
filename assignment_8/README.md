# Assignment 8

The assignment instructions are contained in `Assignment8.pdf`.

## Installation

To compile the Fortran programs, make sure you have `gfortran` installed. The programs were compiled using version 11.2.0.

Python is also required to run the notebooks. The required version is in `.python-version`. If you have `pyenv`
installed, it'll point to that version automatically. To create a virtual environment and install the requirements, run:

```
python -m venv env
source env/bin/activate
pip install -r requirements.txt
```

## Compilation

To compile the modules and programs, run:

```
./compile.sh
```

This stores the programs in the `compiled` directory. To clean out that directory, run:

```
./clean.sh
```

If you can't execute these shell scripts, be sure to make them executable for your user, for example:

```
chmod u+x compile.sh
```

## Execution

The program computes the ground state energy of the 1D Ising model using the real-space renormalization group algorithm. 
To run it, type:

```
compiled/rsrg_ising [--N <N> --max_iter <max_iter> --lambda <lambda> --thres <thres> --diag_method <diag_method> -d/--debug]
```

where the command-line arguments are:

1. `N`: Number of spin sites
2. `max_iter`: Maximum number of iterations to run
3. `lambda`: Coupling between neighboring sites
4. `thres`: Threshold for convergence in terms of successive differences in normalized ground state energy
5. `diag_method`: Method to use to diagonalize (i.e., dsyevr or zheev)
6. `debug`: Whether to display debug information

To run the density matrix renornalization group algorithm, run:

```
compiled/dmrg_ising [--N <N> --max_iter <max_iter> --lambda <lambda> --thres <thres> -d/--debug]
```

where the command-line arguments are interpreted similarly.
