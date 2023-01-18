#include "FFT_2D.hpp"
#include <iostream>
#include <cmath>

class IMAGE_COMPRESSION: protected FFT_2D
{
    public:

        IMAGE_COMPRESSION(){}

        void
        load_image(const char* file_path);

        void
        image_compression(double compression);

        void
        output_image();


    protected:
       
        double compression;

        double compression_factor;

        int
        zeros(cMatrix m);

        void
        quantization(double compression);

        void
        dequantization();


};


