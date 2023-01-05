#include "FFT_2D.hpp"
#include <mpi.h>

int
main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    FFT_2D fft;
    if(mpi_rank == 0)
    fft.generate_random_input(18);
    fft.parallel_solve();
    MPI_Finalize();
    if(mpi_rank == 0){
        fft.evaluate_time_and_error();
    }


    return 0;

}