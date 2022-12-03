#include <iostream>
#include <cmath>
#include <complex>
#include <vector>


using cVector = std::vector<std::complex<double>>;
using Complex = std::complex<double>;
constexpr double pi = 3.14159265358979323846;

cVector recursive_fft ( Complex *p_x, int N, int s)
{
    cVector y;
    y.resize(N);

    if (N==1){
        return {*p_x};
    }

    cVector even;
    cVector odd;

    even = recursive_fft((p_x), N/2, 2*s);
    for(size_t i=0; i<N/2; ++i){
        y[i] = even[i];
    }
    odd = recursive_fft((p_x + s), N/2, 2*s);

    for(size_t i=0; i<N/2; ++i){
        y[i + N/2] = odd[i];
    }

 

    for (size_t k= 0; k < N/2; ++k)
    {
        double kd = static_cast<double> (k);
        double Nd = static_cast<double> (N);
        printf("Before FFT Dimension %d, value at %d: %f and value at %d: %f\n", N,k,  y[k].real(), k+N/2, y[k+N/2].real());
        Complex j(0.0,1.0);
        Complex o = std::exp(-j*2.0*pi*kd/Nd) * y[k + N/2]; 
        Complex p = y[k];
        y[k] = p + o;
        y[k + N/2] = p - o;
        printf("After FFT Dimension %d, value at %d: %f and value at %d: %f\n", N,k,  y[k].real(), k+N/2, y[k+N/2].real());

    }

    return y;

}


int main()
{
    cVector x={{1.0,3.0},{2.0,0.0},{1.5,6.0},{3.0,-1.0},{5.0,2.5},{0.0,0.0},{-1.5,0.0},{7.0,3.5}};

    const int N = x.size();
   
    x = new_fft(&x[0], N,1);

    for(std::size_t i=0; i<N; i++)
        std::cout << x[i].real() << " " << x[i].imag() << std::endl;

    return 0;
}