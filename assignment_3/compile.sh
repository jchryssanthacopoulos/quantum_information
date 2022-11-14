# Compile fortran programs
#
# Example usage:
#   ./compile.sh
#

SRC_DIR="src"
DEST_DIR="compiled"

EXTRA_ARGS="-fno-range-check -J compiled -llapack"


##########################
### Compile Exercise 2 ###
##########################
echo "Compiling Exercise 2 ..."
gfortran $SRC_DIR/exercise_2.f90 -o $DEST_DIR/exercise_2 $EXTRA_ARGS
