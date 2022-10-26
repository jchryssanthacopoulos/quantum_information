# Compile fortran programs
#
# Example usage:
#   ./compile.sh
#

SRC_DIR="src"
DEST_DIR="compiled"

EXTRA_ARGS="-fno-range-check -J compiled"


##########################
### Compile Exercise 1 ###
##########################
echo "Compiling Exercise 1 ..."
gfortran $SRC_DIR/exercise_1.f90 -o $DEST_DIR/exercise_1 $EXTRA_ARGS


##########################
### Compile Exercise 2 ###
##########################
echo "Compiling Exercise 2 ..."
gfortran $SRC_DIR/exercise_2.f90 -o $DEST_DIR/exercise_2 $EXTRA_ARGS


##########################
### Compile Exercise 3 ###
##########################
echo "Compiling Exercise 3 ..."

# compile without optimization
echo "--> With O0 ..."
gfortran $SRC_DIR/exercise_3.f90 -o $DEST_DIR/exercise_3_O0 -O0 $EXTRA_ARGS

# compile with O1 optimization
echo "--> With O1 ..."
gfortran $SRC_DIR/exercise_3.f90 -o $DEST_DIR/exercise_3_O1 -O1 $EXTRA_ARGS

# compile with O2
echo "--> With O2 ..."
gfortran $SRC_DIR/exercise_3.f90 -o $DEST_DIR/exercise_3_O2 -O2 $EXTRA_ARGS

# compile with O3
echo "--> With O3 ..."
gfortran $SRC_DIR/exercise_3.f90 -o $DEST_DIR/exercise_3_O3 -O3 $EXTRA_ARGS

# compile with Os
echo "--> With Os ..."
gfortran $SRC_DIR/exercise_3.f90 -o $DEST_DIR/exercise_3_Os -Os $EXTRA_ARGS

# compile with Ofast
echo "--> With Ofast ..."
gfortran $SRC_DIR/exercise_3.f90 -o $DEST_DIR/exercise_3_Ofast -Ofast $EXTRA_ARGS
