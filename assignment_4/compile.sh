# Compile fortran programs
#
# Example usage:
#   ./compile.sh
#

SRC_DIR="src"
DEST_DIR="compiled"

EXTRA_ARGS="-fno-range-check -J compiled -llapack"


########################
### Compile programs ###
########################
echo "Compiling eigen schrodinger ..."
gfortran $SRC_DIR/eigen_schrodinger.f90 -o $DEST_DIR/eigen_schrodinger $EXTRA_ARGS
