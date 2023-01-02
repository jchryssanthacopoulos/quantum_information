# Compile fortran programs
#
# Example usage:
#   ./compile.sh
#

SRC_DIR="src"
DEST_DIR="compiled"

EXTRA_ARGS="-fno-range-check -J compiled -llapack"


#######################
### Compile Modules ###
#######################
echo "Compiling modules ..."
gfortran -c $SRC_DIR/arg_parse.f90 -o $DEST_DIR/arg_parse.o $EXTRA_ARGS
gfortran -c $SRC_DIR/ising_hamiltonian.f90 -o $DEST_DIR/ising_hamiltonian.o $EXTRA_ARGS


########################
### Compile programs ###
########################
echo "Compiling real-space renormalization group code ..."
gfortran $DEST_DIR/arg_parse.o $DEST_DIR/ising_hamiltonian.o $SRC_DIR/rsrg_ising.f90 -o $DEST_DIR/rsrg_ising $EXTRA_ARGS
