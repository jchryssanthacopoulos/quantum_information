# Compile fortran programs
#
# Example usage:
#   ./compile.sh
#

SRC_DIR="src"
DEST_DIR="compiled"

EXTRA_ARGS="-L/usr/local/lib/lib -fno-range-check -J compiled -llapack -lfftw3"


#######################
### Compile Modules ###
#######################
echo "Compiling modules ..."
gfortran -c $SRC_DIR/arg_parse.f90 -o $DEST_DIR/arg_parse.o $EXTRA_ARGS


########################
### Compile programs ###
########################
echo "Compiling solve time-dependent harmonic oscillator ..."
gfortran $DEST_DIR/arg_parse.o $SRC_DIR/solve_time_dep_ho.f90 -o $DEST_DIR/solve_time_dep_ho $EXTRA_ARGS
