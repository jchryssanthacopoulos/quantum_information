# Assignment 3

The assignment instructions are contained in `Assignment3.pdf`.

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

### Exercise 1

The program multiples two matrices with given dimensions using the specified method. To run the program, type:

```
compiled/exercise_1 [--mat_mul_method <method> --num_rows <nr> --num_cols <nc> --num_inner_dim <ninner> -d/--debug]
```

The optional arguments are

1. `mat_mul_method`: Method to use to multiply the matrices (options are `matmul`, `row-col`, and `col-row`) (default = `matmul`)
2. `num_rows`: Number of rows of the first matrix (default = 10)
3. `num_cols`: Number of columns of the second matrix (default = 10)
4. `num_inner_dim`: Number of inner dimensions (i.e., number of columns of first matrix, rows of second matrix) (default = 10)
5. `debug`: Flag indicating whether to run in debug mode, which prints the inputs and outputs to the screen (default = `False`)

The Python notebook `src/exercise_1.ipynb` runs the program for different matrix sizes and methods, performs fits,
and plots the results.

### Exercise 2

The program computes the eigenvalues of a matrix and saves a histogram of the normalized eigenvalue spacings to a file.
To run it, type:

```
compiled/exericse_2 [--mat_type <mat_type> --output_filename <filename> --ndim <ndim> --nsamples <nsamp> --nbins <nbins> --min_val <min_val> --max_val <max_val> -d/--debug]
```

where the arguments are

1. `mat_type`: Type of matrix to diagonalize (options are `hermitian` and `diag`)
2. `output_filename`: Name of file to save the output histogram
3. `ndim`: Number of rows and columns of the matrix
4. `nsamples`: Number of random matrix to sample to produce histogram
5. `nbins`: Number of bins to use to produce histogram
6. `min_val`: Minimum spacings value of histogram
7. `max_val`: Maximum spacings value of histogram
8. `-d/debug`: Flag indicating whether to run in debug mode, which prints the matrix, eigenvalues, spacings, and their average for each sample

### Exercise 3

The notebook `src/exercise_3.ipynb` can be used to run the program in Exercise 2 and produce fits to the histogram to compare to theory.
