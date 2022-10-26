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

You should see the square root of 5 printed to the screen.

### Exercise 2

To run the program, type:

```
compiled/exercise_2
```

### Exercise 3

To run the program, type:

```
compiled/exercise_3_[opt_flag]
```
where `[opt_flag]` can be `O0`, `O1`, `02`, `03`, `0s`, or `Ofast`.

To plot the run times for the different matrix multiplication methods as a function of the input matrix size, open
the Jupyter notebook `src/plot_exercise_3.ipynb` in your IDE, or start a Jupyter server with

```
jupyter notebook src/plot_exercise_3.ipynb
```
