# Solve Schrodinger equation using different parameters
#
# Usage:
#   ./solve_schrodinger.sh
#


for N in 100 1000 5000; do
    for xmax in 2.5 5 10 15; do
        echo "Running for N = $N, xmin = -$xmax, xmax = $xmax ..."
        compiled/eigen_schrodinger \
            --output_filename data/solution_$N\_$xmax.txt \
            --npoints $N \
            --xmin -$xmax \
            --xmax $xmax
    done
done

echo "Success!"
