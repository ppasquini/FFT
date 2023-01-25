#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
#include "FFT_traits.hpp"

/*!
* Class to perform and analyse the Fast Fourier Trasmorm in 1 dimension
*/
class FFT_1D
{
    public:

        /*!
        * Constructor
        */
        FFT_1D(){}

        /*!
        * Load input from a text file
        * @param file path of the input file. Input vector should be a power of two in order for the programm to work well
        */
        void 
        load_input_from_file(std::string const &file);

        /*!
        * Generate random input of a given dimension
        * @param power dimension of the input will be 2 elevated by the given power
        */
        void
        generate_random_input(unsigned int power);


        /*!
        * Computes the fast fourier trasform in a recursive way
        * @param x input on which compute the fft
        * @param N dimension of the vector
        * @return the result of the fft
        */
        cVector
        recursive_solve(cVector const &x, int N);

        /*!
        * Computes the fast fourier trasform in a iterative way on the input and saves the result in iterative_solution variable
        */
        void 
        iterative_solve();

        /*!
        * Computes the fast fourier trasform in a parallel way on the input using MPI, and saves the result in iterative_solution variable
        */
        void
        parallel_solve();

        /*!
        * Computes the inverse of the fast fourier trasform in a iterative way.
        * @param x vector on which performs the inverse fft
        * @return the result of the computation
        */
        cVector
        inverse_solve(cVector const &vector);

        /*!
        * Calls the test function and then output the result of the parallel solver on the display
        */
        void
        output_and_test();

        /*!
        * Performs the iterative FFT and displays the time taken by both the parallel and the iterative algorithm. Moreover it performs the inverse fft on the solution 
        * of the parallel solver and computes the maximum error with the initial input
        */
        void
        evaluate_time_and_error();

        /*!
        * Saves the result of the parallel solver on a file
        */
        void
        save_output_in_file();

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
        * Computes the discrete fourier trasform on the input and saves the result in discrete_solution
        */
        void
        discrete_solve();        

        /*!
        * Performs the discrete fourier trasform and compares the result with the parallel solution
        */
        void
        test();
        
            
        /*!
        * Dimension of the input
        */
        unsigned int N;

        /*!
        * Time taken by the parallel algorithm
        */
        double time_parallel;

        /*!
        * Time taken by the iterative algorithm
        */
        double time_serial;

        /*!
        * input vector
        */        
        cVector input;

        /*!
        * Solution of the discrete solver
        */
        cVector discrete_solution;

        /*!
        * solution of the parallel solver
        */
        cVector parallel_solution;

        /*!
        * Solution of the iterative solver
        */
        cVector iterative_solution;
};
