#include "FFT_2D.hpp"
#include <iostream>
#include <cmath>

class IMAGE_COMPRESSION: protected FFT_2D
{
    public:

        IMAGE_COMPRESSION(){}

        void
        load_image(const char* file_path);

        double
        image_compression(double compression);

        void
        image_decompression(double comp_facto);

        void
        output_image();

        void
        loadCompression(std::string file_matrix_compressed);

    protected:
       
        double compression;

        double compression_factor;

        cSparseMatrix matrix_compressed;

        void
        quantization(double compression);

        void
        dequantization();


};


