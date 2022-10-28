# Assignment 1

The assignment instructions are contained in `Assignment1.pdf`.

## Installation

To compile the Fortran programs, make sure you have `gfortran` installed. The programs were compiled using version 11.2.0.

Python is also required. The required version is in `.python-version`. If you have `pyenv` installed, it'll point to
that version automatically. To create a virtual environment and install the requirements, run:

```
python -m venv env
source env/bin/activate
pip install -r requirements.txt
```

## Compilation

To compile the programs, run:

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

To run the program, type:

```
compiled/exercise_1
```

This is a test program. You should see the following output:

```
The square root of 5.0 is 2.2361
```

### Exercise 2

To run the program, type:

```
compiled/exercise_2
```

It tests the limits of integer and double precision. You should see the following:

```
The sum of -31616 and 1 using INTEGER*2 is -31615
The sum of 2000000 and 1 using INTEGER*4 is 2000001
The sum of 3.14159278E+32 and 1.41421360E+21 using REAL*4 is 3.14159278E+32
The sum of 3.1415926535897933E+32 and 1.4142135623730950E+21 using REAL*8 is 3.1415926536039354E+32
```

### Exercise 3

To run the program, type:

```
compiled/exercise_3_[opt_flag]
```
where `[opt_flag]` can be `O0`, `O1`, `02`, `03`, `0s`, or `Ofast`.

This program multiplies two matrices using different methods. Enter the size of the input matrices, and the run times
are returned. For example:

```
$ compiled/exercise_3_O0
 Enter matrix dimensions of first matrix:
10, 10
 Enter matrix dimensions of second matrix:
10, 10
Elapsed time for row-col-inner = 1.800000000E-05
Elapsed time for inner-col-row = 1.000000000E-05
Elapsed time for matmul = 3.520000000E-04
```

To plot the run times for the different matrix multiplication methods as a function of the input matrix size, open
the Jupyter notebook `src/plot_exercise_3.ipynb` in your IDE, or start a Jupyter server with

```
jupyter notebook src/plot_exercise_3.ipynb
```

Here's an example plot of run time versus matrix size without optimization:

![image](assignment_1/example_run_times.png}
