#include "FFT_2D.hpp"
#include <chrono>
#include <omp.h>

using namespace std::chrono;

void
FFT_2D::generate_random_input(unsigned int power){
    std::cout << "Loading input..." << std::endl;

    N = std::pow(2,power);
    Nd = static_cast<double>(N);
    input = cMatrix::Random(N,N);
    //Create random test matrix

    std::cout << "Loading done" << std::endl;
    
}

void
FFT_2D::load_input_from_file(std::string file_path){
    std::cout << "Loading input..." << std::endl;

    cSparseMatrix sparseInput;
    Eigen::loadMarket(sparseInput, file_path);
    input = cMatrix(sparseInput);
    N = input.rows();
    Nd = static_cast<double>(N);
    
    std::cout << "Loading done" << std::endl;
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

    std::cout << "Computing iterative FFT" << std::endl;

    cVector input_vector;
    input_vector.resize(N);

    iterative_solution.resize(N,N);

    for (std::size_t i = 0; i < N; i++)
    {
        input_vector = input.row(i);
        iterative_solution.row(i) = iterative_solve_wrapped(input_vector);
    }

    for (std::size_t i = 0; i < N; i++)
    {
        input_vector = iterative_solution.col(i);
        iterative_solution.col(i) = iterative_solve_wrapped(input_vector);
    }
     
}

cVector
FFT_2D::iterative_solve_wrapped(cVector x){
    Complex wd, w, o, p;
    Complex im = {0.0, 1.0};


    //Compute bit reversal
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
    return  x;

}

void
FFT_2D::parallel_solve(){

    const auto t0 = high_resolution_clock::now();

    cVector input_vector;
    input_vector.resize(N);

    parallel_solution.resize(N,N);

    std::cout << "Computing parallel FFT on rows" << std::endl;

    #pragma omp parallel for shared(input, parallel_solution) firstprivate(input_vector) num_threads(num_threads)
    for (std::size_t i = 0; i < N; i++)
    {
        input_vector = input.row(i);
        parallel_solution.row(i) = iterative_solve_wrapped(input_vector);
    }
    
    std::cout << "Computing parallel FFT on columns" << std::endl;

    #pragma omp parallel for shared(parallel_solution) firstprivate(input_vector) num_threads(num_threads)
    for (std::size_t i = 0; i < N; i++)
    {
        input_vector = parallel_solution.col(i);
        parallel_solution.col(i) = iterative_solve_wrapped(input_vector);
    }
    
    const auto t1 = high_resolution_clock::now();
    time_parallel =  duration_cast<milliseconds>(t1 - t0).count();

    std::cout << "Computation complete" << std::endl;
        
}

void
FFT_2D::inverse_fft(){

    cVector input_vector;

    inverse_solution.resize(N, N);

    std::cout << "Computing parallel FFT on rows" << std::endl;

    #pragma omp parallel for shared(inverse_solution, parallel_solution) firstprivate(input_vector) num_threads(num_threads)
    for (std::size_t i = 0; i < N; i++)
    {
        input_vector = parallel_solution.row(i);
        inverse_solution.row(i) = inverse_solve(input_vector);
    }
    
    inverse_solution = (1/(Nd)) * inverse_solution;

    std::cout << "Computing parallel FFT on columns" << std::endl;

    #pragma omp parallel for shared(inverse_solution) firstprivate(input_vector) num_threads(num_threads)
    for (std::size_t i = 0; i < N; i++)
    {
        input_vector = inverse_solution.col(i);
        inverse_solution.col(i) = inverse_solve(input_vector);
    }

    inverse_solution = (1/(Nd)) * inverse_solution;

    std::cout << "Computation complete" << std::endl;

}

cVector
FFT_2D::inverse_solve(cVector x){

    Complex wd, w, o, p;
    Complex im = {0.0, 1.0};

    x = vector_reversal(x, N);

    unsigned int steps = std::log2(N);

    // Iterate Steps
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
    
    return  x;

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

    iterative_solution.resize(0, 0);

    std::cout << "=================================" << std::endl;
    std::cout << "Error Evaluation" <<std::endl << "Inverse fft of the parallel solution compared with the initial input " << std::endl;


    inverse_fft();
    double max_error = ((inverse_solution - input).cwiseAbs()).maxCoeff();

    std::cout << "Max error among all the elements: " << max_error << std::endl;

    inverse_solution.resize(0,0);

}

void
FFT_2D::save_output_in_file(std::string name_file_output){
    cSparseMatrix sparseMatrix = inverse_solution.sparseView();
    Eigen::saveMarket(sparseMatrix, name_file_output);
    sparseMatrix = input.sparseView();
    Eigen::saveMarket(sparseMatrix, "name_file_output");
}


