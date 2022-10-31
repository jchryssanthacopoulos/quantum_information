# Assignment 2

The assignment instructions are contained in `Assignment2.pdf`.

## Installation

To compile the Fortran programs, make sure you have `gfortran` installed. The programs were compiled using version 11.2.0.

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
compiled/exercise_1 [-v/--verbose]
```

The program experiments with debug checkpoints. If the `-v` or `--verbose` option is not provided, the program does
nothing, exiting with

```
 Nothing to print. Exiting ...
```

If the debug flag is provided, a checkpoint is entered and the values of certain variables are printed:

```
 Entering checkpoint ...
Value of variable a is 10
Value of variable b is 20
```

### Exercise 3

To run the program, type:

```
compiled/exercise_3
```

It asks for dimensions of a complex matrix, which it will create using a custom data type called `cmatrix`. It will
initialize it with random data, calculate the trace, and compute the adjoint, displaying all the data to the screen.
It will also write the matrix data into the files `mat.txt` and `mat_adj.txt` (for adjoint) in the current directory.
Here's an example output:

```
 The original matrix is:
.3637 +.7254i   +.4593 +.8790i
.0378 +.9221i   +.5062 +.9254i
The trace of M is .8699 +1.6507i
 Saving to file mat.txt ...
 The adjoint matrix is:
.3637 -.7254i   +.0378 -.9221i
.4593 -.8790i   +.5062 -.9254i
The trace of M is .8699 -1.6507i
 Saving to file mat_adj.txt ...
```
