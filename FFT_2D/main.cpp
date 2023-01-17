#include "FFT_2D.hpp"
#include <mpi.h>

int
main(int argc, char * argv[])
{

    FFT_2D fft;
    fft.load_image();
    fft.image_compression(0.30);
    // fft.evaluate_time_and_error();
    

    return 0;

}