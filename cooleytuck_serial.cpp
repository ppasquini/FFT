#include <iostream>
#include <cmath>
//#include <complex>
//#include <vector>
#include "FFT_traits.hpp"
#include"./FFT_tools.hpp"

//using cVector = std::vector<std::complex<double>>;
//using Complex = std::complex<double>;
//constexpr double pi = 3.14159265358979323846;

int main ()
{
    cVector x;

    const int N = 8;
    double Nd = static_cast<double> (N);

    for (int t = 0; t < N; t++)
    {
        double real = (std::rand() % 10) / 10.0;
        double imag = (std::rand() % 10) / 10.0;
        Complex number = {real, imag};
        x.emplace_back(number);
    }

    Complex im = {0.0,1.0};
    Complex w = std::exp(2.0*pi*im/Nd);

    auto steps = std::log2(N);

    x = FFT::vector_reversal(x,N);

    // Iterate Steps
    for (int i = 1; i <= steps; i++)
    {
        auto power = std::pow(2,i-1);
        // Iterate inside even/odd vectors
        for (int j = 0; j < power; j++)
        {
            if (i != steps)
            {
                for (int k = 0; k < N/(i+1); k++)
                {
                    int s = std::pow(2,i);
                    double jd = static_cast<double> (j);
                    Complex o = std::pow(w,jd) * x[j+k*s+power];
                    x[j+k*s] = x[j+k*s] + o;
                    x[j+k*s+power] = x[j+k*s] - o;
                }
            }
            else
            {
                Complex o = std::pow(w,j) * x[j+power];
                x[j] = x[j] + o;
                x[j+power] = x[j] - o;
            }
        }
    }

    std::cout << "FFT: " << std::endl;
    for (std::size_t i=0; i<N; i++)
    {
        std::cout << x[i].real() << " " << x[i].imag() << std::endl;
    }
}
