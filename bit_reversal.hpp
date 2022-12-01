#include <iostream>
#include <iostream>
#include <vector>
#include <complex>
#include <bits/stdc++.h>


using cVector = std::vector<std::complex<double>>;
using Complex = std::complex<double>;


namespace FFT
{
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
}