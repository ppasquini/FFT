#include "IMAGE_COMPRESSION.hpp"
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_PNG
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

void
IMAGE_COMPRESSION::load_image(const char* file_path){
    std::cout << "Loading image" << std::endl;

    int x,y,n;
    x = 1024;
    y = 1024;
    n = 8;
    unsigned char *data = stbi_load(file_path, &x, &y, &n, 1);

    N = x;

    input.resize(x,y);

    Complex num;

    for(std::size_t i=0; i<y; i++){
        for(std::size_t j=0; j<x; j++){
            num = {static_cast<double>(data[y*i+j]), 0.0};
            input(i,j) = num;
        }
    }

    std::cout << "Done loading" << std::endl;
}

void
IMAGE_COMPRESSION::image_compression(double compression){

    parallel_solve();

    int zeros_before_compression = zeros(parallel_solution);  
   
    quantization(compression);

    int zeros_after_compression = zeros(parallel_solution);    

    std::cout << "Zeros entrys input: " << zeros_before_compression << std::endl;
    std::cout << "Zeros entrys solution: " << zeros_after_compression << " with a percentage of: " << static_cast<double>(zeros_after_compression)/(N*N) * 100 << "%" << std::endl;

    dequantization();

    inverse_fft();

    parallel_solution = 1/(Nd * Nd) * parallel_solution;

    output_image();

}

void
IMAGE_COMPRESSION::quantization(double compression){

    compression = compression / (1.0 - 0.8 * compression);
    double sum = 0.0;
    for(std::size_t i=0; i<N; i++){
        for(std::size_t j=0; j<N; j++){
            sum += std::abs((parallel_solution(i,j)));
        }
    }

    compression_factor = ((sum/(N*N))*(compression));

    for(std::size_t i=0; i<N; i++){ 
        for(std::size_t j=0; j<N; j++){
            parallel_solution(i,j) = (parallel_solution(i,j) / compression_factor);
            if(abs(parallel_solution(i,j)) < 1.0) parallel_solution(i, j) = {0.0, 0.0};  
        }
    }

}

void
IMAGE_COMPRESSION::dequantization(){

    for(std::size_t i=0; i<N; i++){
        for(std::size_t j=0; j<N; j++){
            parallel_solution(i,j) = parallel_solution(i,j) * compression_factor;
        }
    }
}

void
IMAGE_COMPRESSION::output_image(){

    double max, coeff;
    char* v;
    v = (char*) malloc(N*N*sizeof(char));
    
    max = 0;
    for(std::size_t i=0; i<N; i++){
        for(std::size_t j=0; j<N; j++){
            if (std::abs(parallel_solution(i,j)) > max) 
                max = std::abs(parallel_solution(i,j));
        }
    }

    coeff = 255.0 / max;
    
    for(std::size_t i=0; i<N; i++){
        for(std::size_t j=0; j<N; j++){
            v[N*i+j] = static_cast<char>(coeff * std::abs(parallel_solution(i,j)));
        }
    }

    stbi_write_png("output.png", N, N, 1, v, 0);
    free(v);
}

int
IMAGE_COMPRESSION::zeros(cMatrix m){
    int zeros = 0;
    for(std::size_t i=0; i<N; i++){
        for(std::size_t j=0; j<N; j++){
            if(parallel_solution(i, j).real() == 0.0) zeros++;
        }
    }    
    return zeros;
}