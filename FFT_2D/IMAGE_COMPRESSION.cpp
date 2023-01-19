#include "IMAGE_COMPRESSION.hpp"
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_PNG
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

void
IMAGE_COMPRESSION::load_image(const char* file_path, const int size){
    std::cout << "Loading image..." << std::endl;

    int x,y,n;
    x = size;
    y = size;
    n = 8;
    data = stbi_load(file_path, &x, &y, &n, 1);

    N = size;
    Nd = static_cast<double>(N);

    input.resize(x,y);

    Complex num;

    for(std::size_t i=0; i<y; i++){
        for(std::size_t j=0; j<x; j++){
            num = {static_cast<double>(data[y*i+j]), 0.0};
            input.coeffRef(i,j) = num;
        }
    }

    std::cout << "Loading done" << std::endl;
}

void
IMAGE_COMPRESSION::load_image_rgb(const char* file_path, const int size, int channel){

    std::cout << "Loading the ";
    switch (channel)
    {
    case 0:
        std::cout << "red";
        break;
    case 1:
        std::cout << "green";
        break;
    case 2:
        std::cout << "blue";
        break;
    default:
        break;
    }
    std::cout << " channel..." << std::endl;

    int x, y;
    x = size;
    y = size;
    N = size;
    Nd = static_cast<double>(N);

    if(channel == 0){
        int n = 8;
        data = stbi_load(file_path, &x, &y, &n, 3);
        input.resize(x,y);
    }

    Complex num;

    for(std::size_t i=0; i<y; i++){
        for(std::size_t j=0; j<x; j++){
            num = {static_cast<double>(data[3*y*i+j*3+channel]), 0.0};
            input.coeffRef(i,j) = num;
        }
    }

    std::cout << "Loading done" << std::endl;

    if(channel == 2)
        free(data);

}

std::vector<double>
IMAGE_COMPRESSION::image_compression_rgb(const char* file_path, const int size, double compression){

    int zeros_before_compression = 0;
    int zeros_after_compression = 0;

    std::cout << "================================" << std::endl;

    for(size_t color = 0; color < 3; color++){

        std::string scolor = color==0?"r":color==1?"g":"b";

        std::string file_compressed = "Matrix/Matrix_compressed_" + scolor;
        std::string file_not_compressed = "Matrix/Matrix_not_compressed_" + scolor;
       
        Eigen::saveMarket(input, file_not_compressed);

        load_image_rgb(file_path, size, color);

        zeros_before_compression = parallel_solution.size() - parallel_solution.nonZeros();  

        parallel_solve();

        std::cout << "Starting compression..." << std::endl;
   
        quantization(compression);

        matrix_compressed = parallel_solution.sparseView();

        matrix_compressed.makeCompressed();

        Eigen::saveMarket(matrix_compressed, file_compressed); 

        zeros_after_compression += matrix_compressed.size () - matrix_compressed.nonZeros(); 

        compression_factor_rgb.push_back(compression_factor); 

        std::cout << "Compression done: matrix compressed saved in: " << file_compressed << std::endl;

        std::cout << std::endl;

    }
    parallel_solution.resize(0,0);

    std::cout << "================================" << std::endl;

    int memory_saved = (zeros_after_compression - zeros_before_compression) * sizeof(Complex) / 1000; 

    std::cout << "Zeros entrys input: " << zeros_before_compression << std::endl;
    std::cout << "Zeros entrys solution: " << zeros_after_compression << " or the " << static_cast<double>(zeros_after_compression)/(3*N*N) * 100 << "% of the total number of elements" << std::endl;
    std::cout << "Memory saved: " << memory_saved << " KBs" << std::endl;
    std::cout << std::endl;
    
    return compression_factor_rgb;
    
}

double
IMAGE_COMPRESSION::image_compression(double compression){

    Eigen::saveMarket(input, "Matrix/Matrix_not_compressed_bk");

    parallel_solve();

    int zeros_before_compression = parallel_solution.size() - parallel_solution.nonZeros();  

    std::cout << "Starting compression..." << std::endl;
   
    quantization(compression);

    matrix_compressed = parallel_solution.sparseView();
    matrix_compressed.makeCompressed();
    Eigen::saveMarket(matrix_compressed, "Matrix/Matrix_compressed_bk"); 
    parallel_solution.resize(0,0);

    std::cout << "Compression done: matrix compressed saved in: Matrix_compressed" << std::endl;

    std::cout << "================================" << std::endl;

    int zeros_after_compression = matrix_compressed.size () - matrix_compressed.nonZeros();  

    int memory_saved = (zeros_after_compression - zeros_before_compression) * sizeof(Complex) / 1000; 

    std::cout << "Zeros entrys input: " << zeros_before_compression << std::endl;
    std::cout << "Zeros entrys solution: " << zeros_after_compression << " or the " << static_cast<double>(zeros_after_compression)/(3*N*N) * 100 << "% of the total number of elements" << std::endl;
    std::cout << "Memory saved: " << memory_saved << " KBs" << std::endl;
    
    return compression_factor;  
}

void
IMAGE_COMPRESSION::image_decompression(double comp_factor){

    std::cout << "================================" << std::endl;

    std::cout << "Starting decompression..." << std::endl;

    compression_factor = comp_factor;

    dequantization();

    inverse_fft();

    std::cout << "Decompression done: image saved in output.png" << std::endl;

    output_image();
}

void
IMAGE_COMPRESSION::image_decompression_rgb(std::vector<double> comp_factor_rgb, std::vector<std::string> files_matrix_compressed){

    compression_factor_rgb = comp_factor_rgb;

    std::cout << "================================" << std::endl;

    std::cout << "Starting decompression..." << std::endl;

    for(size_t color = 0; color < 3; color ++){

        compression_factor = compression_factor_rgb[color];

        load_compression(files_matrix_compressed[color]);

        dequantization();

        inverse_fft();

        output_image_rgb(color);
    }

    std::cout << "Decompression done: image saved in output.png" << std::endl;
}

void
IMAGE_COMPRESSION::quantization(double compression){

    compression = compression / (1.0 - 0.8 * compression);
    double sum = (parallel_solution.cwiseAbs()).sum();

    compression_factor = ((sum/(N*N))*(compression));

    for(std::size_t i=0; i<N; i++){ 
        for(std::size_t j=0; j<N; j++){
            parallel_solution.coeffRef(i,j) = (parallel_solution.coeff(i,j) / compression_factor);
            if(abs(parallel_solution.coeff(i,j)) < 1.0) parallel_solution.coeffRef(i, j) = {0.0, 0.0};
        }
    }

}

void
IMAGE_COMPRESSION::dequantization(){

    parallel_solution = parallel_solution * compression_factor;
}

void
IMAGE_COMPRESSION::output_image(){

    double max, coeff;
    v = (char*) malloc(N*N*sizeof(char));
    
    max = (inverse_solution.cwiseAbs()).maxCoeff();
    
    coeff = 255.0 / max;
    
    for(std::size_t i=0; i<N; i++){
        for(std::size_t j=0; j<N; j++){
            v[N*i+j] = static_cast<char>(coeff * std::abs(inverse_solution.coeff(i,j)));
        }
    }

    stbi_write_png("output.png", N, N, 1, v, 0);
    free(v);
}

void
IMAGE_COMPRESSION::output_image_rgb(int channel){

    if(channel == 0){
        v = (char*) malloc(3*N*N*sizeof(char));
    }

    
    double max, coeff;
    
    max = (inverse_solution.cwiseAbs()).maxCoeff();

    coeff = 255.0 / max;
    
    
    for(std::size_t i=0; i<N; i++){
        for(std::size_t j=0; j<N; j++){
            v[3*N*i+3*j+channel] = static_cast<char>(coeff * std::abs(inverse_solution.coeff(i,j)));
        }
    }

    if(channel == 2){
        stbi_write_png("output.png", N, N, 3, v, 0);
        free(v);
    }
}

void
IMAGE_COMPRESSION::load_compression(std::string file_matrix_compressed){
    Eigen::loadMarket(matrix_compressed, file_matrix_compressed);
    matrix_compressed.uncompress();
    parallel_solution = cMatrix(matrix_compressed);
    N = matrix_compressed.rows();
    Nd = static_cast<double>(N);

}
