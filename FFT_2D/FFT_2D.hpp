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
        generate_random_input(unsigned int power);

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
        evaluate_time_and_error();

    protected:

        unsigned int 
        bit_reversal(unsigned int value, unsigned int dim);

        cVector
        vector_reversal(cVector x, unsigned int dim);
        
        unsigned int N;

        double Nd = static_cast<double>(N);

        double time_parallel;

        double time_serial;
        
        cMatrix input;

        cVector temp_input;

        cMatrix parallel_solution;

        cMatrix iterative_solution;

        cVector temp_solution;


};
