# Compile fortran programs
#
# Example usage:
#   ./compile.sh
#

SRC_DIR="src"
DEST_DIR="compiled"

EXTRA_ARGS="-fno-range-check -J compiled -llapack"


###########################
##### Compile Modules #####
###########################
echo "Compiling modules ..."
gfortran -c $SRC_DIR/arg_validate.f90 -o $DEST_DIR/arg_validate.o $EXTRA_ARGS
gfortran -c $SRC_DIR/histogram.f90 -o $DEST_DIR/histogram.o $EXTRA_ARGS
gfortran -c $SRC_DIR/mat_ops.f90 -o $DEST_DIR/mat_ops.o $EXTRA_ARGS


##########################
### Compile Exercise 1 ###
##########################
echo "Compiling Exercise 1 ..."
gfortran $DEST_DIR/mat_ops.o $DEST_DIR/arg_validate.o $SRC_DIR/exercise_1.f90 -o $DEST_DIR/exercise_1_O0 $EXTRA_ARGS -O0
gfortran $DEST_DIR/mat_ops.o $DEST_DIR/arg_validate.o $SRC_DIR/exercise_1.f90 -o $DEST_DIR/exercise_1_O1 $EXTRA_ARGS -O1
gfortran $DEST_DIR/mat_ops.o $DEST_DIR/arg_validate.o $SRC_DIR/exercise_1.f90 -o $DEST_DIR/exercise_1_O2 $EXTRA_ARGS -O2
gfortran $DEST_DIR/mat_ops.o $DEST_DIR/arg_validate.o $SRC_DIR/exercise_1.f90 -o $DEST_DIR/exercise_1_O3 $EXTRA_ARGS -O3
gfortran $DEST_DIR/mat_ops.o $DEST_DIR/arg_validate.o $SRC_DIR/exercise_1.f90 -o $DEST_DIR/exercise_1_Ofast $EXTRA_ARGS -Ofast


##########################
### Compile Exercise 2 ###
##########################
echo "Compiling Exercise 2 ..."
gfortran $DEST_DIR/histogram.o $DEST_DIR/mat_ops.o $SRC_DIR/exercise_2.f90 -o $DEST_DIR/exercise_2 $EXTRA_ARGS
