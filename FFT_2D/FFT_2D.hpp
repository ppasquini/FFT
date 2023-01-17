#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include "FFT_traits.hpp"

class FFT_2D
{
    public:

        FFT_2D(){}

        void 
        load_input();

        void
        load_image();

        void
        generate_random_input(unsigned int power);

        void
        image_compression(double compression);

        void 
        iterative_solve();

        void
        iterative_solve_wrapped();

        cVector
        iterative_solve_wrapped(cVector x);

        void
        parallel_solve();

        void
        parallel_solve_wrapped();

        void
        inverse_fft();

        cVector
        inverse_solve(cVector x);

        void
        evaluate_time_and_error();

        void
        output_image();


    protected:

        unsigned int 
        bit_reversal(unsigned int value, unsigned int dim);

        cVector
        vector_reversal(cVector x, unsigned int dim);
        
        unsigned int N;

        double Nd = static_cast<double>(N);

        double time_parallel;

        double time_serial;

        double compression;

        double compression_factor;
        
        cMatrix input;

        cVector temp_input;

        cMatrix parallel_solution;

        cMatrix iterative_solution;

        cVector temp_solution;

        int
        nnz(cMatrix m);

        void
        quantization(double compression);

        void
        dequantization();


};
