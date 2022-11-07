# Compile fortran programs
#
# Example usage:
#   ./compile.sh
#

SRC_DIR="src"
DEST_DIR="compiled"

EXTRA_ARGS="-fno-range-check -J compiled"


###########################
##### Compile Modules #####
###########################
echo "Compiling modules ..."
gfortran -c $SRC_DIR/arg_parse.f90 -o $DEST_DIR/arg_parse.o $EXTRA_ARGS
gfortran -c $SRC_DIR/mat_ops.f90 -o $DEST_DIR/mat_ops.o $EXTRA_ARGS
gfortran -c $SRC_DIR/cmatrix_type.f90 -o $DEST_DIR/cmatrix_type.o $EXTRA_ARGS


##########################
### Compile Exercise 1 ###
##########################
echo "Compiling Exercise 1 ..."
gfortran $DEST_DIR/arg_parse.o $SRC_DIR/exercise_1.f90 -o $DEST_DIR/exercise_1 $EXTRA_ARGS


##########################
### Compile Exercise 2 ###
##########################
echo "Compiling Exercise 2 ..."
gfortran $DEST_DIR/mat_ops.o $DEST_DIR/arg_parse.o $SRC_DIR/exercise_2.f90 -o $DEST_DIR/exercise_2 $EXTRA_ARGS


##########################
### Compile Exercise 3 ###
##########################
echo "Compiling Exercise 3 ..."
gfortran $DEST_DIR/arg_parse.o $DEST_DIR/cmatrix_type.o $SRC_DIR/exercise_3.f90 -o $DEST_DIR/exercise_3 $EXTRA_ARGS
