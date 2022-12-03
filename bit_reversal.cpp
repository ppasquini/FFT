#include <iostream>
#include <vector>
#include <complex>
#include "bit_reversal.hpp"


using cVector = std::vector<std::complex<double>>;
using Complex = std::complex<double>;

 unsigned int bit_reversal( unsigned int value, unsigned int dim){
     unsigned int reversed = 0;

    for(size_t i=0; i<dim;i++){
        reversed <<= 1;
        if((value & 1) == 1) reversed ^= 1;
        value >>= 1;
    }

    return reversed;
}

cVector vector_reversal(cVector x, unsigned int dim){
    cVector y;
    y.resize(x.size());
     unsigned int j = 0;
    for( unsigned int i=0; i<dim; ++i){
        j = bit_reversal(i,std::log2(dim));
        y[j] = x[i];

    }
    return y;

}


int main(){
    cVector vector = {{0.0,7.0},{1.0,6.0},{2.0,5.0},{3.0,4.0},{4.0,3.0},{5.0,2.0},{6.0,1.0},{7.0,0.0}};

    vector = FFT::vector_reversal(vector, vector.size());

    for(std::size_t i=0; i<vector.size(); i++)
        std::cout << vector[i].real() << " " << vector[i].imag() << std::endl;

    

}