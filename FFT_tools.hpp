#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <bits/stdc++.h>


using cVector = std::vector<std::complex<double>>;
using Complex = std::complex<double>;

namespace FFT
{
    constexpr double pi = 3.14159265358979323846;


    cVector recursive_fft ( Complex *p_x, int N, int s)
    {
        cVector y;
        y.resize(N);

        if (N==1){
            return {*p_x};
        }

        cVector even;
        cVector odd;

        even = recursive_fft((p_x), N/2, 2*s);
        for(size_t i=0; i<N/2; ++i){
            y[i] = even[i];
        }
        odd = recursive_fft((p_x + s), N/2, 2*s);

        for(size_t i=0; i<N/2; ++i){
            y[i + N/2] = odd[i];
        }

    

        for (size_t k= 0; k < N/2; ++k)
        {
            double kd = static_cast<double> (k);
            double Nd = static_cast<double> (N);
            Complex j(0.0,1.0);
            Complex o = std::exp(-j*2.0*pi*kd/Nd) * y[k + N/2]; 
            Complex p = y[k];
            y[k] = p + o;
            y[k + N/2] = p - o;

        }

        return y;
    }   

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