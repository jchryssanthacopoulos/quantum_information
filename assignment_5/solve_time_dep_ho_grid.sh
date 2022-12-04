# Solve time-dependent quantum harmonic oscillator with different parameters
#
# Usage:
#   ./solve_time_dep_ho_grid.sh
#


for T in 100 125 150 175 200 225 250 275 300 325 350 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500; do
    echo "Running for T = $T ..."
    compiled/solve_time_dep_ho \
        --output_filename data/psi_Nt_1000_tmax_$T.txt \
        --num_t_pts 1000 \
        --tmax $T
done

echo "Success!"
