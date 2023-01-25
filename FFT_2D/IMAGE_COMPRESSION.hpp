#include "FFT_2D.hpp"
#include <iostream>
#include <cmath>

/*!
* Class create to compress both black and white and colored images using the 2 dimensional fft implemented in the FFT_2D class
*/
class IMAGE_COMPRESSION: protected FFT_2D
{
    public:

        /*!
        * Constructor of the class
        * @param threads number of threads used for OPENMP
        */
        IMAGE_COMPRESSION(int threads):FFT_2D(threads) {}

        /*!
        * Load image from a file. For black and white images
        * @param file path of the input file.
        * @param size dimension of the image. It should be a power of two to work well
        */
        void
        load_image(std::string const &file_path, const int size);

        /*!
        * Function to compress black and white images. It saves in the Matrix folder the matrix before and after the compression.
        * @param compression parameter given for the compression
        * @return the compression factor calculated using the mean of the matrix
        */
        double
        image_compression(double compression);

        /*!
        * Function to compress colored images. It saves in the Matrix folder the matrix before and after the compression for every color channel
        * @param file_path file path of the image to be compressed
        * @param size size of the image to be compressed. It should be a power of two to work well
        * @param compression parameter given for the compression
        * @return vector of compression factors calculated using the mean of the matrix for every color channel
        */
        std::vector<double>
        image_compression_rgb(std::string const &file_path, const int size, double compression);

        /*!
        * Function to decompress black and white images
        * @param comp_facto compression factor returned by the compression function
        * @return the compression factor calculated using the mean of the matrix
        */
        void
        image_decompression(double comp_facto);

        /*!
        * Function to decompress colored images
        * @param comp_factor_rgb vector of the compression factors returned by the compression function for every color channel
        * @param files_matrix_compressed vector containing the files of the compressed matrixes for every color channel
        */
        void
        image_decompression_rgb(std::vector<double> const &comp_factor_rgb, std::vector<std::string> const &files_matrix_compressed);
        
        /*!
        * Loads a compressed matrix from a file
        * @param file_matrix_compressed path of the file with the matrix compressed
        */
        void
        load_compression(std::string const &file_matrix_compressed);          
 
    protected:

        /*!
        * Function to divide the matix obtained by the FFT by the compression factor. It helps to get a sparse matrix
        * @param compression desired
        */
        void
        quantization(double compression);

        /*!
        * Function to multiply the matix obtained by the FFT by the given factor. It helps recreate the original image from the sparse one
        */
        void
        dequantization();

        /*!
        * Load image from a file. For colored images
        * @param file path of the input file.
        * @param size dimension of the image. It should be a power of two to work well
        */
        void
        load_image_rgb(std::string const &file_path, const int size, int channel);

        /*!
        * Function to save the decompressed image
        */
        void
        output_image();

        /*!
        * Function to save the decompressed image
        * @param channel color channel computed at the moment
        */
        void
        output_image_rgb(int channel);      
       
        /*!
        * Compression desired
        */      
        double compression;

        /*!
        * Factor calculated and used for the compression
        */
        double compression_factor;

        /*!
        * Vector of the compression factors calculated for every color channel
        */
        std::vector<double> compression_factor_rgb;

        /*!
        * Pointer to the variable used to store the image read from input
        */
        std::unique_ptr<unsigned char*> input_data;


        /*!
        * Pointer to the variable used to store a matrix in order to write it as an image
        */
        std::unique_ptr<std::string> output;

        /*!
        * Compressed matrix saved as a sparse matrix
        */
        cSparseMatrix matrix_compressed;
};


