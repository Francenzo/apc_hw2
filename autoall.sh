
rm *.stdout

make
sbatch auto-bridges-serial
sbatch auto-bridges-pthreads16
sbatch auto-bridges-openmp16
sbatch auto-bridges-mpi16
