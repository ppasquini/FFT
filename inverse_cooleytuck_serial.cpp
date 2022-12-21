#include <iostream>
#include <cmath>
#include <mpi.h>
#include <chrono>
#include "FFT_traits.hpp"
#include "./FFT_tools.hpp"

cVector inverse_cooleytuck_serial(cVector x){

    unsigned int N = x.size();
    Complex wd, w, o, p;
    Complex im = {0.0, 1.0};

    //Compute bit reversal
    x = FFT::vector_reversal(x, N);

    unsigned int steps = std::log2(N);

    // Iterate Steps
    for (int i = 1; i <= steps; i++)
    {
        auto power = std::pow(2,i);
        // Calculate primitive root w
        wd = std::exp(2.0*FFT::pi*im/power);
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

    for(int index = 0; index < N; ++index){
        x[index] /= N;
    }

    return x;
}

cVector cooleytuck_serial(cVector x){

    unsigned int N = x.size();
    Complex wd, w, o, p;
    Complex im = {0.0, 1.0};

    //Compute bit reversal
    x = FFT::vector_reversal(x, N);

    unsigned int steps = std::log2(N);

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

    return x;
}


int main (int argc, char* argv[])
{


    cVector x;
    cVector iterative_solution;
    cVector inverse_solution;
    //Number of elements
    const unsigned int N = std::pow(2,18);
    x.resize(N);


    //Create random test vector
    for (int t = 0; t < N; t++)
    {
        double real = (std::rand() % 100) / 20.0 - 2.5;
        double imag = (std::rand() % 100) / 20.0 - 2.5;

        x[t] =  {real, imag};
    }





    //Find the result of the DTF using standard algorithm
    iterative_solution = cooleytuck_serial(x);

    inverse_solution = inverse_cooleytuck_serial(iterative_solution);

    std::cout << "IFFT: " << std::endl;
    
    bool correct = true;
    
    //Print results
    for (std::size_t i=0; i<N; i++)
    {
        std::cout << "At index " << i << ": " << std::endl;
        std::cout << "Right solution: " << x[i].real() << " " << x[i].imag() << std::endl;
        std::cout << "Inverse solution: " << inverse_solution[i].real() << " " << inverse_solution[i].imag() << std::endl;
        std::cout << "Iterative solution: " << iterative_solution[i].real() << " " << iterative_solution[i].imag() << std::endl << std::endl;
    }

    //Evaluate correctness by comparing with the result of the recursive FFT
    double tol = 1.e-6;
    std::cout.precision(16);
    if(inverse_solution.size() != x.size()) correct = false;
    for (std::size_t i=0; i<N; i++)
    {         

        if(!((std::abs(inverse_solution[i].real() - x[i].real()) < tol) && 
                (std::abs(inverse_solution[i].imag() - x[i].imag()) < tol))){
            std::cout << "Value wrong: " << inverse_solution[i] << " at index: " << i << ". Right value: " << x[i] << std::fixed << std::endl << std::endl;
            correct = false;
        }

    }    
    if(correct) 
        std::cout << "Algorithm completed successfully" << std::endl;
    else 
        std::cout << "Something is wrong!" << std::endl;

    return 0;
}
