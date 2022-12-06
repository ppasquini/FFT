#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

#include"./FFT_tools.hpp"

using cVector = std::vector<std::complex<double>>;
using Complex = std::complex<double>;

int main ()
{
    cVector x;
    //Number of elements
    const int N = 8;
    double Nd = static_cast<double> (8);
    //w coefficient
    Complex wd;

    //Create random vector
    for (int t = 0; t < N; t++)
    {
        double real = (std::rand() % 10) / 10.0;
        double imag = (std::rand() % 10) / 10.0;
        Complex number = {real, imag};
        x.emplace_back(number);
    }

    Complex im = {0.0,1.0};

    //Number of steps of the algorithm
    auto steps = std::log2(N);

    //Find the result of the FFT using the recursive algorithm
    cVector recursive_solution = FFT::recursive_fft(x, N);

    //COmpute vit reversal
    x = FFT::vector_reversal(x,N);

    // Iterate Steps
    for (int i = 1; i <= steps; i++)
    {

        auto power = std::pow(2,i);
        //Calculate coefficient w
        wd = std::exp(2.0*FFT::pi*im/power);
        // Iterate inside even/odd vectors
        for (int j = 0; j < power/2; j++)
        {
            for (int k = j; k < N; k+=power)
            {
                double jd = static_cast<double> (j);
                Complex o = std::pow(wd, jd) * x[k + power/2];
                Complex p = x[k]; 
                x[k] = p + o;
                x[k + power/2] = p - o;
            }
        }
    }

    std::cout << "FFT: " << std::endl;

    //Print results
    for (std::size_t i=0; i<N; i++)
    {
        std::cout << x[i].real() << " " << x[i].imag() << std::endl;
        std::cout << recursive_solution[i].real() << " " << recursive_solution[i].imag() << std::endl << std::endl;

    }
    
    //Evaluate correctness by comparing with the result of the recursive FFT
    bool correct = true;
    double tol = 1.e-14;
    if(x.size() != recursive_solution.size()) correct = false;
    for (std::size_t i=0; i<N; i++)
    {         
        size_t l;                                      
        for(l=0; l<N; ++l){
            if((std::abs(x[i].real() - recursive_solution[l].real()) < tol) && (std::abs(x[i].imag() - recursive_solution[l].imag())<tol)) break;
        }
        if(l == N){
            std::cout << "Value wrong: " << x[i] << std::endl;
            correct = false;
        }

    }    
    if(correct) std::cout << "Both algorithm match!" << std::endl;
    else std::cout << "Something is wrong!" << std::endl;

}