#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

#include "./FFT_tools.hpp"

using cVector = std::vector<std::complex<double>>;
using Complex = std::complex<double>;

int main ()
{
    cVector x;
    //Number of elements
    const int N = std::pow(2,20);
    x.resize(N);
    //w roots and temporary variables
    Complex wd, w, o, p;


    //Create random test vector
    for (int t = 0; t < N; t++)
    {
        double real = (std::rand() % 100) / 20.0 - 2.5;
        x[t] =  {real, 1.0};
    }

    Complex im = {0.0,1.0};

    //Number of steps of the algorithm
    auto steps = std::log2(N);

    //Find the result of the FFT using the recursive algorithm
    cVector recursive_solution = FFT::recursive_fft(x, N);

    //Find the result of the DTF using standard algorithm
    cVector discrete_solution = FFT::dft(x, N);

    //Compute bit reversal
    x = FFT::vector_reversal(x,N);

    // Iterate Steps
    for (int i = 1; i <= steps; i++)
    {
        auto power = std::pow(2,i);
        // Calculate primitive root w
        wd = std::exp(-2.0*FFT::pi*im/power);
        w = {1.0,0.0};
        // Iterate inside even/odd vectors
        for (int j = 0; j < power/2; j++) // j = exponent of w
        {
            for (int k = j; k < N; k+=power)
            {
                o = w * x[k + power/2];
                p = x[k]; 
                x[k] = p + o;
                x[k + power/2] = p - o;
            }
            w *= wd;
        }
    }

    std::cout << "FFT: " << std::endl;

    //Print results
    for (std::size_t i=0; i<N; i++)
    {
        std::cout << "Non-recursive solution: " << x[i].real() << " " << x[i].imag() << std::endl;
        std::cout << "Recursive solution: " << recursive_solution[i].real() << " " << recursive_solution[i].imag() << std::endl;
        std::cout << "Discrete solution: " << discrete_solution[i].real() << " " << discrete_solution[i].imag() << std::endl << std::endl;
    }
    
    //Evaluate correctness by comparing with the result of the recursive FFT
    bool correct = true;
    double tol = 1.e-10;
    if(x.size() != discrete_solution.size()) correct = false;
    for (std::size_t i=0; i<N; i++)
    {         

        if(!((std::abs(x[i].real() - discrete_solution[i].real()) < tol) && 
                (std::abs(x[i].imag() - discrete_solution[i].imag()) < tol))){
            std::cout << "Value wrong: " << x[i] << " Index: " << i << std::endl;
            correct = false;
        }

    }    
    if(correct) 
        std::cout << "Both algorithm match!" << std::endl;
    else 
        std::cout << "Something is wrong!" << std::endl;

    return 0;
}