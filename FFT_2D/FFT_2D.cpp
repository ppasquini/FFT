#include "FFT_2D.hpp"
#include <mpi.h>
#include <chrono>
#include <omp.h>
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_PNG
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


using namespace std::chrono;

void
FFT_2D::generate_random_input(unsigned int power){
    std::cout << "Loading input" << std::endl;

    N = std::pow(2,power);
    input = cMatrix::Random(N,N);
    //Create random test matrix

    std::cout << "Done loading" << std::endl;
    
}

void
FFT_2D::load_image(){
    std::cout << "Loading image" << std::endl;

    int x,y,n;
    x = 256;
    y = 256;
    n = 8;
    unsigned char *data = stbi_load("test.png", &x, &y, &n, 1);

    N = x;

    input.resize(x,y);

    Complex num;

    for(std::size_t i=0; i<y; i++){
        for(std::size_t j=0; j<x; j++){
            num = {static_cast<double>(data[y*i+j]), 0.0};
            input(i,j) = num;
        }
    }

    std::cout << "Done loading" << std::endl;
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
        y(j) = x(i);

    }
    return y;
}


void
FFT_2D::iterative_solve(){

    std::cout << "Starting fft iterative" << std::endl;

    cVector input_vector;
    input_vector.resize(N);

    iterative_solution.resize(N,N);


    for (std::size_t i = 0; i < N; i++)
    {
        cVector input_vector = input.row(i);
        iterative_solution.row(i) = iterative_solve_wrapped(input_vector);
    }

    for (std::size_t i = 0; i < N; i++)
    {
        cVector input_vector = input.col(i);
        iterative_solution.col(i) = iterative_solve_wrapped(input_vector);
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


cVector
FFT_2D::iterative_solve_wrapped(cVector x){
    Complex wd, w, o, p;
    Complex im = {0.0, 1.0};

    //std::cout << "Iterative fft" << std::endl;

    //Compute bit reversal
    //std::cout << "Bit Reversal" << std::endl;
    //std::cout << "Thread: " << omp_get_thread_num() << " Vector size " << x.size() << std::endl;
    x = vector_reversal(x, N);

    unsigned int steps = std::log2(N);

    // Iterate Steps
    //std::cout << "Computation" << std::endl;
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
                o = w * x.coeffRef(k + power/2);
                p = x.coeffRef(k); 
                x.coeffRef(k) = p + o;
                x.coeffRef(k + power/2) = p - o;
            }
            w *= wd;
        }

    }
    //std::cout << "end" << std::endl;
    return  x;

}

void
FFT_2D::parallel_solve(){

    const auto t0 = high_resolution_clock::now();


    cVector input_vector;
    input_vector.resize(N);

    parallel_solution.resize(N,N);

    std::cout << "Starting fft parallel on rows" << std::endl;

    #pragma omp parallel for shared(input, parallel_solution) firstprivate(input_vector) num_threads(6)
    for (std::size_t i = 0; i < N; i++)
    {
        //std::cout << "Thread: " << omp_get_thread_num() << " Row: " << i << std::endl;
        cVector input_vector = input.row(i);
        //std::cout << "Thread: " << omp_get_thread_num() << " End fft on row: " << i << std::endl;
        parallel_solution.row(i) = iterative_solve_wrapped(input_vector);
        //std::cout << "Thread: " << omp_get_thread_num() << " End copy" << std::endl;
    }
    
    std::cout << "Starting fft parallel on columns" << std::endl;

    #pragma omp parallel for shared(parallel_solution) firstprivate(input_vector) num_threads(6)
    for (std::size_t i = 0; i < N; i++)
    {
        //std::cout << "Thread: " << omp_get_thread_num() << " Row: " << i << std::endl;
        cVector input_vector = parallel_solution.col(i);
        //std::cout << "Thread: " << omp_get_thread_num() << " End fft on row: " << i << std::endl;
        parallel_solution.col(i) = iterative_solve_wrapped(input_vector);
        //std::cout << "Thread: " << omp_get_thread_num() << " End copy" << std::endl;
    }
    
    const auto t1 = high_resolution_clock::now();
    time_parallel =  duration_cast<milliseconds>(t1 - t0).count();

    std::cout << "Done computation" << std::endl;
        
}

void
FFT_2D::inverse_fft(){

    cVector input_vector;

    std::cout << "Starting fft parallel on rows" << std::endl;

    #pragma omp parallel for shared(parallel_solution) firstprivate(input_vector) num_threads(6)
    for (std::size_t i = 0; i < N; i++)
    {
        //std::cout << "Thread: " << omp_get_thread_num() << " Row: " << i << std::endl;
        cVector input_vector = parallel_solution.row(i);
        //std::cout << "Thread: " << omp_get_thread_num() << " End fft on row: " << i << std::endl;
        parallel_solution.row(i) = inverse_solve(input_vector);
        //std::cout << "Thread: " << omp_get_thread_num() << " End copy" << std::endl;
    }
    
    std::cout << "Starting fft parallel on columns" << std::endl;

    #pragma omp parallel for shared(parallel_solution) firstprivate(input_vector) num_threads(6)
    for (std::size_t i = 0; i < N; i++)
    {
        //std::cout << "Thread: " << omp_get_thread_num() << " Row: " << i << std::endl;
        cVector input_vector = parallel_solution.col(i);
        //std::cout << "Thread: " << omp_get_thread_num() << " End fft on row: " << i << std::endl;
        parallel_solution.col(i) = inverse_solve(input_vector);
        //std::cout << "Thread: " << omp_get_thread_num() << " End copy" << std::endl;
    }

    std::cout << "Done computation" << std::endl;

}

cVector
FFT_2D::inverse_solve(cVector x){

    Complex wd, w, o, p;
    Complex im = {0.0, 1.0};

    x = vector_reversal(x, N);

    unsigned int steps = std::log2(N);

    // Iterate Steps
    //std::cout << "Computation" << std::endl;
    for (int i = 1; i <= steps; i++)
    {
        auto power = std::pow(2,i);
        // Calculate primitive root w
        wd = std::exp(2.0*pi*im/power);
        w = {1.0,0.0};
        // Iterate inside even/odd vectors
        for (int j = 0; j < power/2; j++) // j = exponent of w
        {
            for (int k = j; k < N; k+=power)
            {
                o = w * x.coeffRef(k + power/2);
                p = x.coeffRef(k); 
                x.coeffRef(k) = p + o;
                x.coeffRef(k + power/2) = p - o;
            }
            w *= wd;
        }

    }
    //std::cout << "end" << std::endl;
    return  x;

}

void
FFT_2D::image_compression(double compression){

    parallel_solve();

    parallel_solution = 1/Nd * parallel_solution;

    inverse_fft();

    parallel_solution = 1/Nd * parallel_solution;

    output_image();

}

void
FFT_2D::output_image(){

    double max, coeff;
    char* v;
    v = (char*) malloc(N*N*sizeof(char));
    
    max = 0;
    for(std::size_t i=0; i<N; i++){
        for(std::size_t j=0; j<N; j++){
            if (std::abs(parallel_solution(i,j)) > max) 
                max = std::abs(parallel_solution(i,j));
        }
    }

    coeff = 255.0 / max;
    
    for(std::size_t i=0; i<N; i++){
        for(std::size_t j=0; j<N; j++){
            v[N*i+j] = static_cast<char>(coeff * std::abs(parallel_solution(i,j)));
        }
    }

    stbi_write_png("output.png", N, N, 1, v, 0);
    free(v);
}

/*
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

*/

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


