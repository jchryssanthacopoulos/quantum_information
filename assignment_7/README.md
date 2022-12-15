# Assignment 7

The assignment instructions are contained in `Assignment7.pdf`.

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

The program generates the Hamiltonian of the 1D Ising model and solves for the eigenvalues and eigenvectors. To run it, 
type:

```
compiled/ising_model [--N <N> --lambda <lambda> --output_filename <output_filename> --debug <debug_level>]
```

where the command-line arguments are:

1. `N`: Number of spin sites
2. `lambda`: Coupling between neighboring sites
3. `output_filename`: Name of file to save solution
4. `debug`: Whether to display debug information
