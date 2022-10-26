# Compile fortran programs
#
# Example usage:
#   ./compile.sh
#


##########################
### Compile Exercise 1 ###
##########################
echo "Compiling Exercise 1 ..."
gfortran src/exercise_1.f90 -J compiled -o compiled/exercise_1


##########################
### Compile Exercise 2 ###
##########################
echo "Compiling Exercise 2 ..."
gfortran src/exercise_2.f90 -J compiled -o compiled/exercise_2


##########################
### Compile Exercise 3 ###
##########################
echo "Compiling Exercise 3 ..."

# compile without optimization
echo "--> With O0 ..."
gfortran src/exercise_3.f90 -J compiled -o compiled/exercise_3_O0 -O0

# compile with O1 optimization
echo "--> With O1 ..."
gfortran src/exercise_3.f90 -J compiled -o compiled/exercise_3_O1 -O1

# compile with O2
echo "--> With O2 ..."
gfortran src/exercise_3.f90 -J compiled -o compiled/exercise_3_O2 -O2

# compile with O3
echo "--> With O3 ..."
gfortran src/exercise_3.f90 -J compiled -o compiled/exercise_3_O3 -O3

# compile with Os
echo "--> With Os ..."
gfortran src/exercise_3.f90 -J compiled -o compiled/exercise_3_Os -Os

# compile with Ofast
echo "--> With Ofast ..."
gfortran src/exercise_3.f90 -J compiled -o compiled/exercise_3_Ofast -Ofast
