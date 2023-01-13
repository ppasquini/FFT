#include "FFT_2D.hpp"
#include <mpi.h>

int
main(int argc, char * argv[])
{

    FFT_2D fft;
    fft.generate_random_input(10);
    fft.parallel_solve();
    fft.evaluate_time_and_error();
    


    return 0;

}