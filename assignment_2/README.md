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
compiled/exercise_1 [-d/--debug]
```

The program experiments with debug checkpoints. If the `-d` or `--debug` option is not provided, the program does
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

### Exercise 2

To run the program, type:

```
compiled/exercise_2 [-d/--debug]
```

This program multiplies two matrices using different methods. Enter the size of the input matrices, and the run times
are returned. For example:

```
$ compiled/exercise_2
 Enter number of rows, columns, and inner dimension:
10, 10, 10
Elapsed time for matmul = 3.3100000000E-04
Max abs error for row-col = 8.8817841970E-16
Elapsed time for row-col = 8.0000000000E-06
Max abs error for col-row = 8.8817841970E-16
Elapsed time for col-row = 7.0000000000E-06
```

If your inputs are not positive integers, you'll get the following errors:

```
# non-integers
$ compiled/exercise_2
 Enter number of rows, columns, and inner dimension:
a b c
 Dimensions need to be integers!

# non-positive integers
$ compiled/exercise_2
 Enter number of rows, columns, and inner dimension:
1, 2, -1
 Dimensions must be greater than zero!
```

When you run with the `-d/--debug` flag, it prints the input and output matrices as an aid to debugging:

```
 Running in debug mode ...
 Matrix A = 
   0.82   0.79
   0.80   0.99
 Matrix B = 
   0.25   0.71
   0.48   0.96
 Product using matmul = 
   0.58   1.35
   0.67   1.53
Elapsed time for matmul = 1.7000000000E-05
 Product using row-col = 
   0.58   1.35
   0.67   1.53
Max abs error for row-col = 0.0000000000E+00
Elapsed time for row-col = 3.0000000000E-06
 Product using col-row = 
   0.58   1.35
   0.67   1.53
Max abs error for col-row = 0.0000000000E+00
Elapsed time for col-row = 2.0000000000E-06
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
+0.6862 +0.6708i   +0.4096 +0.0576i
+0.8869 +0.4427i   +0.1159 +0.7512i
 The trace of M is +0.8021 +1.4219i
 Saving to file mat.txt ...
 The adjoint matrix is:
+0.6862 -0.6708i   +0.8869 -0.4427i
+0.4096 -0.0576i   +0.1159 -0.7512i
 The trace of M is +0.8021 -1.4219i
 Saving to file mat_adj.txt ...
 Deleting matrices ...
```
