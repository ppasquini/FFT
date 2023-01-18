#include "IMAGE_COMPRESSION.hpp"
#include <mpi.h>
#include <chrono>

int
main(int argc, char * argv[])
{

    IMAGE_COMPRESSION compressor;
    if(argc > 2){
        compressor.load_image(argv[1]);
        float compression = atof(argv[2]);
        double compression_factor = compressor.image_compression(compression);
        IMAGE_COMPRESSION new_compressor;
        compressor.loadCompression("Matrix_compressed");
        compressor.image_decompression(compression_factor);
    }
    else{
        std::cout << "Not enough input! Please enter the name of the path of the image and a paramenter with the compression desired" << std::endl;
        return 1;
    }    // fft.evaluate_time_and_error();
    


    return 0;

}