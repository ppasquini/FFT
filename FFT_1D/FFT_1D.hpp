#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include "FFT_traits.hpp"

class FFT_1D
{
    public:

        FFT_1D(){}

        void 
        load_input();

        void
        generate_random_input(unsigned int power);

        void
        discrete_solve();

        cVector
        recursive_solve(cVector x, int N);

        void 
        iterative_solve();

        void
        parallel_solve();

        void
        test();

        void
        output_and_test();

        void
        evaluate_time_and_error();

    protected:

        unsigned int 
        bit_reversal(unsigned int value, unsigned int dim);

        cVector
        vector_reversal(cVector x, unsigned int dim);
        
        unsigned int N;

        double time_parallel;

        double time_serial;
        
        cVector input;

        cVector discrete_solution;

        cVector parallel_solution;

        cVector iterative_solution;




};
