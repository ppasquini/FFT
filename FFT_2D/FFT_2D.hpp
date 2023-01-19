#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include "FFT_traits.hpp"

class FFT_2D
{
    public:

        FFT_2D(int threads):num_threads{threads}{}

        void
        generate_random_input(unsigned int power);

        void
        load_input_from_file(std::string file_path);

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
        save_output_in_file(std::string name_file_output);

    protected:

        unsigned int 
        bit_reversal(unsigned int value, unsigned int dim);

        cVector
        vector_reversal(cVector x, unsigned int dim);

        int
        num_threads;
        
        unsigned int N;

        double Nd;

        double time_parallel;

        double time_serial;
        
        cMatrix input;

        cVector temp_input;

        cMatrix parallel_solution;

        cMatrix iterative_solution;

        cMatrix inverse_solution;

        cVector temp_solution;

};
