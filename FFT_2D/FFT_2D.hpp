#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include "FFT_traits.hpp"

/*!
* Class to perform and analyse the Fast Fourier Trasmorm in 1 dimension
*/
class FFT_2D
{
    public:

        /*!
        * Constructor of the class
        * @param threads number of threads used for OPENMP
        */
        FFT_2D(int threads):num_threads{threads}{}

        /*!
        * Generate random input of a given dimension
        * @param power dimension of the input will be 2 elevated by the given power
        */
        void
        generate_random_input(unsigned int power);

        /*!
        * Load input from a file
        * @param file path of the input file. Input matrix shoud be a power of two in order for the programm to work well
        */
        void
        load_input_from_file(std::string const &file_path);
       
       /*!
        * Computes the fast fourier trasform on every row and then on every column of the input matrix and saves the result in iterative_solution variable
        */
        void 
        iterative_solve();


       /*!
        * Computes the fast fourier trasform in a iterative way on the given vector and returns the result
        * @param x vector on which performs the iterative fft
        * @return result of the computation
        */
        cVector
        iterative_solve_wrapped(cVector const &vector);

       /*!
        * Computes the iterative fast fourier trasform in a parallel way using OPENMP on every row and then on every column of the input matrix and saves the result
        * in iterative_solution variable
        */
        void
        parallel_solve();

       /*!
        * Computes the inverse fast fourier trasform in a parallel way using OPENMP on every row and then on every column of the input matrix and saves the result
        * in inverse_solution variable
        */
        void
        inverse_fft();

       /*!
        * Computes the iterative fast fourier trasform on a given vector and return the result
        * @param x vector on which performs the inverse FFT
        * @return result of the computation
        */
        cVector
        inverse_solve(cVector const &vector);

        /*!
        * Performs the iterative FFT and displays the time taken by both the parallel and the iterative algorithm. Moreover it performs the inverse fft on the solution 
        * of the parallel solver and computes the maximum error with the initial input
        */
        void
        evaluate_time_and_error();

        /*!
        * Saves the result of the parallel solver on a file of a given path
        * @param name_file_output name of the file where it saves the result
        */
        void
        save_output_in_file(std::string const &name_file_output);

    protected:

        /*!
        * Performs the bit reversal algorithm on a single value
        * @param value value on which perform the bit reversal
        * @param dim dimension in log2 of the value (es dimension of 8 will be 3)
        * @return the value with the bit reversed
        */
        unsigned int 
        bit_reversal(unsigned int value, unsigned int dim);

        /*!
        * Performs the bit reversal algorithm on the indices of a given vector
        * @param x vector on which perform the bit reversal
        * @param dim dimension of the vector
        */
        cVector
        vector_reversal(cVector const &x, unsigned int dim);

        /*!
        * Number of threads
        */
        int
        num_threads;

        /*!
        * Dimensio of the matrix
        */        
        unsigned int N;

        /*!
        * Dimension of the matrix in double
        */
        double Nd;

        /*!
        * Time taken by the parallel algorithm
        */
        double time_parallel;

        /*!
        * Time taken by the iterative algorithm
        */
        double time_serial;

        /*!
        * Input matrix
        */       
        cMatrix input;

        /*!
        * solution of the parallel solver
        */
        cMatrix parallel_solution;

        /*!
        * solution of the iterative solver
        */
        cMatrix iterative_solution;

        /*!
        * solution of the inverse FFT solver
        */
        cMatrix inverse_solution;

};
