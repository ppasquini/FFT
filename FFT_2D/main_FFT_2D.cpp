#include "FFT_2D.hpp"
#include <mpi.h>
#include <chrono>

int
main(int argc, char * argv[])
{

    if(argc > 1){
        int threads = atoi(argv[1]);
        FFT_2D fft_2D(threads);
        std::string file(argv[2]);
        if(file =="random" && argc > 2){
             int exp = atoi(argv[3]);
             fft_2D.generate_random_input(exp);
        }
        else{
             fft_2D.load_input_from_file(file);
        }
        fft_2D.parallel_solve();
        fft_2D.evaluate_time_and_error();
        //fft_2D.save_output_in_file("matrix_output");
    }
    else{
        std::cout << "Not enough input! Please enter the name of the file contaning the input or the word random followed by a dimension to create a random input of the given dimension" << std::endl;
        return 1;
    }


    return 0;

}