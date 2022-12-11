# Assignment 6

The assignment instructions are contained in `Assignment6.pdf`.

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

The program generates the density matrix for a many-body quantum system and its left and right reduced density matrices.
To run it,
type:

```
compiled/density_matrix [--N <N> --D <D> --type <type> --M <M> --output_filename <output_filename> --debug <debug_level>]
```

where the command-line arguments are:

1. `N`: Number of subsystems in the quantum state
2. `dim`: Number of dimensions of each subsystem
3. `type`: Type of system to compute (options are "separable", "bell", and "generic")
4. `M`: Number of subsystems in the right system
5. `output_filename`: Name of file to save density matrices
6. `debug`: Debug level (options are 0, 1, and 2)

As an example, to compute the density matrices for a 6-body system in 4 dimensions, type:

```
$ compiled/density_matrix --N 6 --D 4 --M 3 --type generic --output_filename data/density_mat_N6_D4_M3.txt
 N = 6
 D = 4
 type = generic
 M = 3
 debug = 0
 output_filename = data/density_mat_N6_D4_M3.txt                     
Trace of density matrix = 1.0000 +0.0000i
Trace of left reduced density matrix = 1.0000 +0.0000i
Trace of right reduced density matrix = 1.0000 +0.0000i
Entropy of left partition = 3.67859
Entropy of right partition = 3.67859
```
