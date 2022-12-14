#include "FFT_2D.hpp"
#include <mpi.h>
#include <chrono>

using namespace std::chrono;

void
FFT_2D::generate_random_input(unsigned int power){
    std::cout << "Loading input" << std::endl;

    N = std::pow(2,power);
    input.resize(N);
    //Create random test matrix
    for (int k = 0; k < N; k++)
    {
        input[k].resize(N);

        for (int t = 0; t < N; t++)
        {
            double real = (std::rand() % 100) / 20.0 - 2.5;
            double imag = (std::rand() % 100) / 20.0 - 2.5;

            input[k][t] =  {real, imag};
        }
    }
    
}

unsigned int 
FFT_2D::bit_reversal(unsigned int value, unsigned int dim){
    // Number reversed
    unsigned int reversed = 0;

    for(size_t i=0; i<dim;i++){
        //Shift the number to the left
        reversed <<= 1;
        //If value as a 1 for the first bit from the right, write 1 on the first bit of reversed
        if((value & 1) == 1) reversed ^= 1;
        //Shift the value of 1 bit to the right
        value >>= 1;
    }

    return reversed;
}

cVector
FFT_2D::vector_reversal(cVector x, unsigned int dim){
    cVector y;
    y.resize(x.size());
    unsigned int j = 0;
    //Save every elements on its reversed position
    for( unsigned int i=0; i<dim; ++i){
        j = bit_reversal(i,std::log2(dim));
        y[j] = x[i];

    }
    return y;
}

void
FFT_2D::iterative_solve(){

    temp_input.resize(N);

    for (std::size_t i = 0; i < N; i++)
    {
        std::copy(input[i].begin(), input[i].end(), temp_input.begin());
        iterative_solve_wrapped();
        std::copy(temp_solution.begin(), temp_solution.end(), iterative_solution[i].begin());
    }

    for (std::size_t i = 0; i < N; i++)
        for (std::size_t c = 0; c < N; c++)
            iterative_solution[i][c] = 1/Nd * iterative_solution[i][c];


    for (std::size_t i = 0; i < N; i++)
    {
        for (std::size_t c = 0; c < N; c++)
            temp_input[c] = iterative_solution[c][i];
        iterative_solve_wrapped();
        for (std::size_t c = 0; c < N; c++)
            iterative_solution[c][i] = 1/Nd * temp_solution[c];
    }
        
}

void
FFT_2D::iterative_solve_wrapped(){
    Complex wd, w, o, p;
    Complex im = {0.0, 1.0};

    //Compute bit reversal
    temp_solution = vector_reversal(temp_input, N);

    unsigned int steps = std::log2(N);

    // Iterate Steps
    for (int i = 1; i <= steps; i++)
    {
        auto power = std::pow(2,i);
        // Calculate primitive root w
        wd = std::exp(-2.0*pi*im/power);
        w = {1.0,0.0};
        // Iterate inside even/odd vectors
        for (int j = 0; j < power/2; j++) // j = exponent of w
        {
            for (int k = j; k < N; k+=power)
            {
                o = w * temp_solution[k + power/2];
                p = temp_solution[k]; 
                temp_solution[k] = p + o;
                temp_solution[k + power/2] = p - o;
            }
            w *= wd;
        }
    }
}

void
FFT_2D::parallel_solve(){

    const auto t0 = high_resolution_clock::now();

    temp_input.resize(N);

    for (std::size_t i = 0; i < N; i++)
    {
        std::copy(input[i].begin(), input[i].end(), temp_input.begin());
        parallel_solve_wrapped();
        std::copy(temp_solution.begin(), temp_solution.end(), iterative_solution[i].begin());
    }

    for (std::size_t i = 0; i < N; i++)
    {
        for (std::size_t c = 0; c < N; c++)
            temp_input[c] = parallel_solution[c][i];
        parallel_solve_wrapped();
        for (std::size_t c = 0; c < N; c++)
            parallel_solution[c][i] = 1/Nd * temp_solution[c];
    }

    const auto t1 = high_resolution_clock::now();
    time_parallel =  duration_cast<milliseconds>(t1 - t0).count();
        
}

void
FFT_2D::parallel_solve_wrapped(){

    int mpi_rank,mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);


    Complex wd, w, o, p;
    const Complex im = {0.0,1.0};

    if(mpi_rank == 0){
        std::cout << "=================================" << std::endl;
        std::cout << "Solving prolem parallel" << std::endl;
        std::cout << "PHASE: I" << std::endl;
        //PHASE I: Compute bit reversal
        temp_solution = vector_reversal(temp_input, N);
    }



    unsigned int N = temp_solution.size();
    MPI_Bcast(&N, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    const unsigned int local_size = N / mpi_size;

    cVector local_vector;
    local_vector.resize(local_size);

    MPI_Scatter(temp_solution.data(), local_size , MPI_DOUBLE_COMPLEX, local_vector.data(), local_size, MPI_DOUBLE_COMPLEX , 0, MPI_COMM_WORLD);

    std::cout << "Thread: " << mpi_rank << " in PHASE: II" << std::endl;

    //PHASE II: compute with local elements
    unsigned int initial_steps = std::log2(N) - std::log2(mpi_size);
    for(size_t i = 1; i <= initial_steps; i++){

        auto power = std::pow(2,i);
        // Calculate primitive root w
        wd = std::exp(-2.0*pi*im/power);
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
        wd = std::exp(-2.0*pi*im/power);
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
        //std::cout << "Thread: " << mpi_rank << " send to " << swap << std::endl;
        MPI_Sendrecv(local_vector.data(), local_size, MPI_DOUBLE_COMPLEX, swap, 0, swap_vector.data(), local_size, MPI_DOUBLE_COMPLEX, swap, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  

        //std::cout << "Thread: " << mpi_rank << " Step: " << i << std::endl;

        // Iterate inside even/odd vectors
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

    MPI_Gather(local_vector.data(), local_size, MPI_DOUBLE_COMPLEX,
                temp_solution.data(), local_size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    std::cout << "Thread: " << mpi_rank << " DONE" << std::endl;
}

void
FFT_2D::evaluate_time_and_error(){
    std::cout << "=================================" << std::endl;
    std::cout << "Time Evaluation" <<std::endl << "Compared with best serial implementation(iterative): " << std::endl;   

    const auto t0 = high_resolution_clock::now();
    iterative_solve();
    const auto t1 = high_resolution_clock::now();
    time_serial = duration_cast<milliseconds>(t1 - t0).count();

    std::cout << "Time taken by Serial Implementation: " << time_serial << " ms" << std::endl << "Time taken by Parallel Implementation: " << time_parallel << " ms" << std::endl;
    auto difference = time_serial - time_parallel;
    std::cout << "Time gained: "<< difference << " ms" << std::endl;

    std::cout << "=================================" << std::endl;
    std::cout << "Error Evaluation" <<std::endl << "Compared with best serial implementation(iterative): " << std::endl;

    //TODO: COMPLETE ERROR COMPUTATION
    //double relative_error = std::abs(parallel_solution. - iterative_solution)/std::abs(iterative_solution);

    //std::cout << "Relative error: " << relative_error << std::endl;
}


