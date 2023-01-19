#include "IMAGE_COMPRESSION.hpp"
#include <mpi.h>
#include <chrono>

int
main(int argc, char * argv[])
{


    if(argc > 5){
        int threads = atoi(argv[1]);
        IMAGE_COMPRESSION compressor(threads);

        const char* file_path = argv[2];
        int dimension_image = atoi(argv[3]);
        float compression = atof(argv[4]);
        std::string color = argv[5];

        if(color == "bk"){
            compressor.load_image(file_path, dimension_image);
            double compression_factor = compressor.image_compression(compression);
            compressor.load_compression("Matrix/Matrix_compressed_bk");
            compressor.image_decompression(compression_factor);
        }
        else{
            std::vector<double> compression_factors; 
            compression_factors = compressor.image_compression_rgb(file_path, dimension_image, compression);
            std::vector<std::string> compression_files = {"Matrix/Matrix_compressed_r", "Matrix/Matrix_compressed_g", "Matrix/Matrix_compressed_b"};
            compressor.image_decompression_rgb(compression_factors, compression_files);
        }
    }
    else{
        std::cout << "Not enough inputs! Please enter:\n- the number of threads\n- name of the image\n- the size of its side\n- paramenter from 0 to 1 for the desired compression\n- \"nk\" for a black and white image or \"rgb\" for a colored one" << std::endl;
        return 1;
    }   
    


    return 0;

}