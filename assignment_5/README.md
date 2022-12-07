# Assignment 5

The assignment instructions are contained in `Assignment5.pdf`.

## Installation

To compile the Fortran programs, make sure you have `gfortran` installed. The programs were compiled using version 11.2.0.

You also need to install the [FFTW library](https://www.fftw.org/). For example, on Unix machines, you can download the source
code, then run:

```
./configure --prefix <lib_path>
make
make install
```

where `lib_path` is the directory to install the library.

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

To change the library path where it searches for the FFTW library, modify the `-L` flag in `EXTRA_ARGS` in the compile
script.

## Execution

The program generates the wavefunction of the time-dependent quantum harmonic oscillator. To run it,
type:

```
compiled/solve_time_dep_ho [--xmin <xmin> --xmax <xmax> --tmax <tmax> --num_x_pts <num_x_pts> --num_t_pts <num_t_pts> --output_filename <output_filename> -d/--debug]
```

where the command-line arguments are:

1. `xmin`: Minimum x value in the domain to solve for (default = `-5.0`)
2. `xmax`: Maximum x value (default = `5.0`)
3. `tmax`: Maximum t value to simulate (default = `10.0`)
4. `num_x_pts`: Number of points to discretize x coordinate (default = `1000`)
5. `num_t_pts`: Number of points to discretize t coordinate (default = `100`)
6. `output_filename`: Name of output file to save wavefunction (default = `solution.txt`)
7. `debug`: Whether to display debug information to the screen

As an example:

```
$ compiled/solve_time_dep_ho --xmin -15 --xmax 15 --tmax 10 --num_x_pts 100 --num_t_pts 100
 xmin = -15.000 
 xmax = 15.000  
 tmax = 10.000  
 num_x_pts = 100     
 num_t_pts = 100     
 debug = F       
 output_filename = solution.txt
Elapsed time to evolve state = 2.9180000000E-03
```

To run several discretization schemes in batch, run

```
./solve_time_dep_ho_grid.sh
```

To plot the wavefunction, run the notebook `src/plot_wavefunction.ipynb`.
