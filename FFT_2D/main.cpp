#include "IMAGE_COMPRESSION.hpp"
#include <mpi.h>
#include <chrono>

int
main(int argc, char * argv[])
{

    IMAGE_COMPRESSION compressor;
    IMAGE_COMPRESSION new_compressor;
    if(argc > 3){
        /*
        compressor.load_image(argv[1], atoi(argv[2]));
        float compression = atof(argv[3]);
        double compression_factor = compressor.image_compression(compression);
        compressor.loadCompression("Matrix_compressed.txt");
        compressor.image_decompression(compression_factor);*/ //for black and white

        std::vector<double> compression_factors; 
        compression_factors = compressor.image_compression_rgb(argv[1], atoi(argv[2]), atof(argv[3]));
        std::vector<std::string> compression_files = {"Matrix_compressed_r", "Matrix_compressed_g", "Matrix_compressed_b"};
        std::this_thread::sleep_for(std::chrono::seconds(5));
        compressor.image_decompression_rgb(compression_factors, compression_files);
    }
    else{
        std::cout << "Not enough inputs! Please enter the name of the image, the size of its side and a paramenter from 0 to 1 for the desired compression" << std::endl;
        return 1;
    }    // fft.evaluate_time_and_error();
    


    return 0;

}