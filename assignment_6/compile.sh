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
gfortran -c $SRC_DIR/entropy.f90 -o $DEST_DIR/entropy.o $EXTRA_ARGS


########################
### Compile programs ###
########################
echo "Compiling density matrix code ..."
gfortran $DEST_DIR/arg_parse.o $DEST_DIR/entropy.o $SRC_DIR/density_matrix.f90 -o $DEST_DIR/density_matrix $EXTRA_ARGS
