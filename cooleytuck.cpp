#include <iostream>
#include <cmath>
#include <complex>
#include <vector>


using cVector = std::vector<std::complex<double>>;
using Complex = std::complex<double>;
constexpr double pi = 3.14159265358979323846;

void fft (cVector& x, int N)
{
    if (N==1)
        return;
    
    cVector even;
    cVector odd;

    for (std::size_t i=0; i<N; i++)
    {
        if (i%2==0)
            even.emplace_back(x[i]);
        else
            odd.emplace_back(x[i]);
    }

    fft(even, N/2);
    fft(odd, N/2);

    for (std::size_t k=0; k<N/2; k++)
    {
        double kd = static_cast<double> (k);
        double Nd = static_cast<double> (N);
        Complex j(0.0,1.0);
        Complex o = std::exp(-j*2.0*pi*kd/Nd) * odd[k]; 
        x[k] = even[k] + o;
        x[k+N/2] = even[k] - o;
    }

}

int main()
{
    cVector x={{1.0,0.0},{2.0,0.0},{-1.5,0.0},{3.0,0.0},{-0.5,0.0},{0.0,0.0},{1.0,0.0},{-2.5,0.0}};

    const int N = 8;
   
    fft(x,N);

    for(std::size_t i=0; i<N; i++)
        std::cout << x[i].real() << " " << x[i].imag() << std::endl;

    return 0;
}