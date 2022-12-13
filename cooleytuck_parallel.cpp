#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <mpi.h>

#include "./FFT_tools.hpp"

using cVector = std::vector<std::complex<double>>;
using Complex = std::complex<double>;

cVector cooleytuck_parallel(cVector vector){
    int mpi_rank,mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    Complex wd, w, o, p;
    const Complex im = {0.0,1.0};



    unsigned int N = vector.size();
    MPI_Bcast(&N, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    const unsigned int local_size = N / mpi_size;

    cVector local_vector;
    local_vector.resize(local_size);

    MPI_Scatter(vector.data(), local_size , MPI_DOUBLE_COMPLEX, local_vector.data(), local_size, MPI_DOUBLE_COMPLEX , 0, MPI_COMM_WORLD);

    std::cout << "Thread: " << mpi_rank << " in PHASE: II" << std::endl;

    //PHASE II: compute with local elements
    unsigned int initial_steps = std::log2(N) - std::log2(mpi_size);
    for(size_t i = 1; i <= initial_steps; i++){

        auto power = std::pow(2,i);
        // Calculate primitive root w
        wd = std::exp(-2.0*FFT::pi*im/power);
        w = {1.0,0.0};
        // Iterate inside even/odd vectors
        for (size_t j = 0; j < power/2; j++) // j = exponent of w
        {
            for (size_t k = j; k < local_size; k+=power)
            {
                o = w * local_vector[k + power/2];
                p = local_vector[k]; 
                local_vector[k] = p + o;
                local_vector[k + power/2] = p - o;
            }
            w *= wd;
        }        

    }


    std::cout << "Thread: " << mpi_rank << " in PHASE: III" << std::endl;

    //PHASE III: swap elements and  compute
    unsigned int final_steps = std::log2(N);
    cVector swap_vector;
    swap_vector.resize(local_size);

    for(size_t i = initial_steps + 1; i <= final_steps; i++){

        auto power = std::pow(2,i);
        // Calculate primitive root w
        wd = std::exp(-2.0*FFT::pi*im/power);
        int sign;
        int size = static_cast <unsigned int> ((power/(local_size))); //numero vettori scambio
        // Determine if the processor has to send data "forwards" or "backwards"
        if((mpi_rank % size) < (size/2)){
            sign = 1;
            w = std::pow(wd, (mpi_rank % size) * local_size);
        }
        else{
            sign = -1;
            w = std::pow(wd, ((mpi_rank + sign * (size/2)) % size) * local_size);
        }
        //Rank of processor communicating with
        unsigned int swap = mpi_rank + sign * size/2;
        std::cout << "Thread: " << mpi_rank << " send to " << swap << std::endl;
        MPI_Send(local_vector.data() , local_size, MPI_DOUBLE_COMPLEX, swap, 0, MPI_COMM_WORLD);
        MPI_Recv(swap_vector.data(), local_size , MPI_DOUBLE_COMPLEX, swap, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        std::cout << "Thread: " << mpi_rank << " Step: " << i << std::endl;

        // Iterate inside even/odd vectorsS
        for (size_t k = 0; k < local_size; k++)
            {
            if(sign > 0){
                o = w * swap_vector[k];
                p = local_vector[k]; 
                local_vector[k] = p + o;                    
            }
            else{
                o = w * local_vector[k];
                p = swap_vector[k]; 
                local_vector[k] = p - o;
            }
            w *= wd;
        }      
    }


    cVector result;
    result.resize(N);

    MPI_Gather(local_vector.data(), local_size, MPI_DOUBLE_COMPLEX,
                result.data(), local_size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    std::cout << "Thread: " << mpi_rank << " DONE" << std::endl;


    return result;


}


int main (int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    cVector x;
    cVector recursive_solution;
    cVector discrete_solution;
    //Number of elements
    const unsigned int N = std::pow(2,11);
    x.resize(N);

    if (mpi_rank == 0){
        //Create random test vector
        for (int t = 0; t < N; t++)
        {
            double real = (std::rand() % 100) / 20.0 - 2.5;
            double imag = (std::rand() % 100) / 20.0 - 2.5;

            x[t] =  {real, imag};
        }
    

        //Find the result of the FFT using the recursive algorithm
        recursive_solution = FFT::recursive_fft(x, N);

        //Find the result of the DTF using standard algorithm
        discrete_solution = FFT::dft(x, N);

        

        std::cout << "PHASE: I" << std::endl;
        //PHASE I: Compute bit reversal
        x = FFT::vector_reversal(x,N);

    }

    x = cooleytuck_parallel(x);

    MPI_Barrier(MPI_COMM_WORLD);

    if(mpi_rank == 0){
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
        std::cout.precision(16);
        if(x.size() != discrete_solution.size()) correct = false;
        for (std::size_t i=0; i<N; i++)
        {         

            if(!((std::abs(x[i].real() - discrete_solution[i].real()) < tol) && 
                    (std::abs(x[i].imag() - discrete_solution[i].imag()) < tol))){
                std::cout << "Value wrong: " << x[i] << " at index: " << i << ". Discrete value: " << discrete_solution[i] << ". Recursive value: " << recursive_solution[i] << std::fixed << std::endl;
                correct = false;
            }

        }    
        if(correct) 
            std::cout << "Both algorithm match!" << std::endl;
        else 
            std::cout << "Something is wrong!" << std::endl;
    }

    MPI_Finalize();
    return 0;
}