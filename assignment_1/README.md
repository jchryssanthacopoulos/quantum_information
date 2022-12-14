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
where `[opt_flag]` can be `O0`, `O1`, `02`, `03`, `0s`, or `Ofast`. More information on optimization flags can be found [here](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html).

This program multiplies two matrices using different methods. Enter the size of the input matrices, and the run times
are returned. For example:

```
$ compiled/exercise_3_O0
 Enter number of rows, columns, and inner dimension:
10, 10, 10
Elapsed time for matmul = 3.3100000000E-04
Max abs error for row-col-inner = 8.8817841970E-16
Elapsed time for row-col-inner = 8.0000000000E-06
Max abs error for inner-col-row = 8.8817841970E-16
Elapsed time for inner-col-row = 7.0000000000E-06
```

To plot the run times for the different matrix multiplication methods as a function of the input matrix size, open
the Jupyter notebook `src/plot_exercise_3.ipynb` in your IDE, or start a Jupyter server with

```
jupyter notebook src/plot_exercise_3.ipynb
```

Here's an example plot of run time versus matrix size without optimization:

![image](example_run_times.png)

## Running on CloudVeneto

To run the code on CloudVeneto, you first need an account. After obtaining one, follow the instructions [here](https://userguide.cloudveneto.it/en/latest/GettingStarted.html#) to set up a keypair, which is a secret key for interacting with a virtual machine after you create it. If you already have an SSH
key, you can import it instead. Then follow the instructions for creating a virtual machine.

You can only log into your virtual machine by SSHing into the gateway. Before doing that, copy your private SSH key
onto the gateway machine. For example, to copy it to the root directory, type:

```
scp /path/to/private/key [username]@gate.cloudveneto.it:~
```

Then SSH into the gateway using your username and password:

```
ssh [username]@gate.cloudveneto.it
```

From there, you can SSH into your VM:

```
ssh -i /path/to/private/key ubuntu@[VM_IP_address]
```

where `[VM_IP_address]` can be found on the web dashboard.
