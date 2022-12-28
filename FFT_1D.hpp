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
        output();


    protected:
        unsigned int 
        bit_reversal(unsigned int value, unsigned int dim);

        cVector
        vector_reversal(cVector x, unsigned int dim);
        
        unsigned int N;

        cVector input;

        cVector discrete_solution;

        cVector solution;


};
