#include "FFT_1D.hpp"
#include <mpi.h>

int
main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    FFT_1D fft;
    if(mpi_rank == 0)
    fft.load_input();
    fft.parallel_solve();
    MPI_Finalize();
    if(mpi_rank == 0)
    fft.output();


    return 0;




}