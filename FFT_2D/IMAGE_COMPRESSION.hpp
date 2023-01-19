#include "FFT_2D.hpp"
#include <iostream>
#include <cmath>

class IMAGE_COMPRESSION: protected FFT_2D
{
    public:

        IMAGE_COMPRESSION(){}

        void
        load_image(const char* file_path, const int size);

        void
        load_image_rgb(const char* file_path, const int size, int channel);

        double
        image_compression(double compression);

        std::vector<double>
        image_compression_rgb(const char* file_path, const int size, double compression);

        void
        image_decompression(double comp_facto);

        void
        image_decompression_rgb(int color);

        void
        output_image();

        void
        output_image_rgb(int channel);

        void
        loadCompression(std::string file_matrix_compressed);

        void
        loadCompression_rgb(std::vector<double> comp_factor_rgb, std::vector<std::string> files_matrix_compressed);
 

    protected:
       
        double compression;

        double compression_factor;

        std::vector<double> compression_factor_rgb;

        unsigned char *data;

        char* v;

        cSparseMatrix matrix_compressed;

        void
        quantization(double compression);

        void
        dequantization();


};


