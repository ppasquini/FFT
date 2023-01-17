#include "FFT_2D.hpp"
#include <mpi.h>

int
main(int argc, char * argv[])
{

    FFT_2D fft;
    if(argc > 1){
        fft.load_image(argv[1]);
        int compression = atoi(argv[2]);
    }
    else{
        std::cout << "Not enough input! Please enter the name of the path of the image and a paramenter with the compression desired" << std::endl;
        return 1;
    }
    fft.image_compression(0.5);
    // fft.evaluate_time_and_error();
    

    return 0;

}