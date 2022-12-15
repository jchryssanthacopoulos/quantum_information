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


########################
### Compile programs ###
########################
echo "Compiling Ising model code ..."
gfortran $DEST_DIR/arg_parse.o $SRC_DIR/ising_model.f90 -o $DEST_DIR/ising_model $EXTRA_ARGS
