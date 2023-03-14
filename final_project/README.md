# Time-evolving Block Decimation with Matrix Product States

This package implements the TEBD algorithm for MPS.

## Installation

The required Python version is in `.python-version`. If you have `pyenv` installed, it'll point to that version
automatically. To create a virtual environment and install the package, run:

```
python -m venv env
source env/bin/activate
pip install -e .
```

## Execution

The time evolution step is given by the following tensor network:

![image](figures/apply_gate_step.png)

This contracted tensor network represents the local density matrix:

![image](figures/local_density_matrix.png)

The average energy of two sites is:

![image](figures/average_energy.png)
