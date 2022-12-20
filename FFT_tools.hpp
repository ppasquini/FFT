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

    cVector dft(cVector x, int dim){
        Complex i = {0.0,1.0};
        cVector result;
        result.resize(dim);

        for(size_t k=0; k<dim; ++k){
            Complex sum = 0;;
            for(size_t j = 0; j < dim; j++){
                sum += x[j] * (cos((2*pi*k*j)/dim) - i * sin((2*pi*k*j)/dim));
            }
            result[k] = sum;
        }
        
        return result;
    }



    cVector recursive_fft ( cVector x, int N)
    {
        // Stopping condition
        if (N==1){
            return x;
        }

        // Preparing necessary vectors and values
        double Nd = static_cast<double> (N);
        Complex j = {0.0,1.0};
        Complex wn = std::exp(-j*2.0*pi/Nd);
        Complex w = {1.0,0.0};

        cVector y;
        y.resize(N);

        cVector y_even;
        cVector y_odd;
        y_even.resize(N/2);
        y_odd.resize(N/2);

        cVector even;
        cVector odd;
        even.resize(N/2);
        odd.resize(N/2);

        // Save even elements in even vector and odd elements in odd vector
        for (std::size_t i = 0; i < N; i++)
        {
            if (i%2)
                odd[i/2] = x[i];
            else
                even[i/2] = x[i];
        }

        y_even = recursive_fft(even, N/2);
        y_odd = recursive_fft(odd, N/2);

    
        for (size_t k = 0; k < N/2; ++k)
        {
            Complex o =  w * y_odd[k]; 
            Complex p = y_even[k];
            y[k] = p + o;
            y[k + N/2] = p - o;
            w *= wn; 
        }

        return y;
    }   

    unsigned int bit_reversal( unsigned int value, unsigned int dim){
        // Number reversed
        unsigned int reversed = 0;

        for(size_t i=0; i<dim;i++){
            //Shift the number to the left
            reversed <<= 1;
            //If value as a 1 for the first bit from the right, write 1 on the first bit of reversed
            if((value & 1) == 1) reversed ^= 1;
            //Shift the value of 1 bit to the right
            value >>= 1;
        }

        return reversed;
    }

    cVector vector_reversal(cVector x, unsigned int dim){
        cVector y;
        y.resize(x.size());
        unsigned int j = 0;
        //Save every elements on its reversed position
        for( unsigned int i=0; i<dim; ++i){
            j = bit_reversal(i,std::log2(dim));
            y[j] = x[i];

        }
        return y;

    }

        cVector parallel_vector_reversal(cVector x, unsigned int dim){
        cVector y;
        y.resize(x.size());
        unsigned int j = 0;
        //Save every elements on its reversed position
        #pragma omp parallel firstprivate(j)
        {
            #pragma omp for 
            for( unsigned int i=0; i<dim; ++i){
                j = bit_reversal(i,std::log2(dim));
                y[j] = x[i];

            }
        }
        return y;

    }
}