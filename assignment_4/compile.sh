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
echo "Compiling eigen schrodinger ..."
gfortran $DEST_DIR/arg_parse.o $SRC_DIR/eigen_schrodinger.f90 -o $DEST_DIR/eigen_schrodinger $EXTRA_ARGS
