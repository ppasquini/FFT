#include "FFT_1D.hpp"
#include <mpi.h>

int
main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    FFT_1D fft;
    if(mpi_rank == 0){
        if(argc > 1){
            std::string file(argv[1]);
            if(file =="random" && argc > 2){
                int exp = atoi(argv[2]);
                fft.generate_random_input(exp);
            }
            else{
                fft.load_input_from_file(file);
            }
        }
        else{
            std::cout << "Not enough input! Please enter the name of the file contaning the input or the word random followed by a dimension to create a random input of the given dimension" << std::endl;
            return 1;
        }
    }
    fft.parallel_solve();
    MPI_Finalize();
    if(mpi_rank == 0){
        //fft.output_and_test(); //Only use for small input
        fft.evaluate_time_and_error();
        fft.save_output_in_file();
    }


    return 0;




}