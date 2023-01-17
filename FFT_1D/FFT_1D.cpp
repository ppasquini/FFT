#include "FFT_1D.hpp"
#include <mpi.h>
#include <chrono>


using namespace std::chrono;

void
FFT_1D::generate_random_input(unsigned int power){
    std::cout << "Generating random input of dimension: " << std::pow(2, power) << std::endl;

    N = std::pow(2,power);
    input.resize(N);
    //Create random test vector
    for (int t = 0; t < N; t++)
    {
        double real = (std::rand() % 100) / 20.0 - 2.5;
        double imag = (std::rand() % 100) / 20.0 - 2.5;

        input[t] =  {real, imag};
    }
}

void
FFT_1D::load_input_from_file(std::string file){
    std::cout << "Loading input from file: " << file << std::endl;
    std::ifstream infile(file);

    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        double real, imag;
        if (!(iss >> real >> imag)) { 
            std::cout << "Error loading" << std::endl;
            break; 
        } // error

        input.push_back({real, imag});
    }

    N = input.size();

    infile.close();
}

void
FFT_1D::discrete_solve(){
    Complex i = {0.0,1.0};
    discrete_solution.resize(N);

    for(size_t k=0; k<N; ++k){
        Complex sum = 0;;
        for(size_t j = 0; j < N; j++){
            sum += input[j] * (cos((2*pi*k*j)/N) - i * sin((2*pi*k*j)/N));
        }
        discrete_solution[k] = sum;
    }
}

unsigned int 
FFT_1D::bit_reversal(unsigned int value, unsigned int dim){
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
FFT_1D::vector_reversal(cVector x, unsigned int dim){
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

cVector
FFT_1D::recursive_solve(cVector x, int N){
    // Stopping condition
    if (N==1){
        return x;
    }

    // Preparing necessary vectors and values
    double Nd = static_cast<double> (N);
    Complex j = {0.0,1.0};
    Complex wn = std::exp(-j*2.0*pi/Nd);
    Complex w = {1.0,0.0};

    cVector y;
    y.resize(N);

    cVector y_even;
    cVector y_odd;
    y_even.resize(N/2);
    y_odd.resize(N/2);

    cVector even;
    cVector odd;
    even.resize(N/2);
    odd.resize(N/2);

    // Save even elements in even vector and odd elements in odd vector
    for (std::size_t i = 0; i < N; i++)
    {
        if (i%2)
            odd[i/2] = x[i];
        else
            even[i/2] = x[i];
    }

    y_even = FFT_1D::recursive_solve(even, N/2);
    y_odd = FFT_1D::recursive_solve(odd, N/2);


    for (size_t k = 0; k < N/2; ++k)
    {
        Complex o =  w * y_odd[k]; 
        Complex p = y_even[k];
        y[k] = p + o;
        y[k + N/2] = p - o;
        w *= wn; 
    }

    return y;
}

void
FFT_1D::iterative_solve(){
    Complex wd, w, o, p;
    Complex im = {0.0, 1.0};

    //Compute bit reversal
    iterative_solution = vector_reversal(input, N);

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
                o = w * iterative_solution[k + power/2];
                p = iterative_solution[k]; 
                iterative_solution[k] = p + o;
                iterative_solution[k + power/2] = p - o;
            }
            w *= wd;
        }
    }
}

void
FFT_1D::parallel_solve(){

    const auto t0 = high_resolution_clock::now();

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
        parallel_solution = vector_reversal(input, N);
    }



    unsigned int N = parallel_solution.size();
    MPI_Bcast(&N, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    const unsigned int local_size = N / mpi_size;

    cVector local_vector;
    local_vector.resize(local_size);

    MPI_Scatter(parallel_solution.data(), local_size , MPI_DOUBLE_COMPLEX, local_vector.data(), local_size, MPI_DOUBLE_COMPLEX , 0, MPI_COMM_WORLD);

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
                parallel_solution.data(), local_size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    if(mpi_rank == 0){
        const auto t1 = high_resolution_clock::now();
        time_parallel =  duration_cast<milliseconds>(t1 - t0).count();
    }

    std::cout << "Thread: " << mpi_rank << " DONE" << std::endl;
}

cVector
FFT_1D::inverse_solve(cVector x){

    unsigned int N = x.size();
    Complex wd, w, o, p;
    Complex im = {0.0, 1.0};

    //Compute bit reversal
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



void
FFT_1D::test(){
    std::cout << "=================================" << std::endl;
    std::cout << "Testing output:" << std::endl;

    bool correct = true;
    //Evaluate correctness by comparing with the result of the recursive FFT
    double tol = 1.e-6;
    std::cout.precision(16);
    if(parallel_solution.size() != discrete_solution.size()) correct = false;
    for (std::size_t i=0; i<N; i++)
    {         

        if(!((std::abs(parallel_solution[i].real() - discrete_solution[i].real()) < tol) && 
                (std::abs(parallel_solution[i].imag() - discrete_solution[i].imag()) < tol))){
            std::cout << "Value wrong: " << parallel_solution[i] << " at index: " << i << ". Discrete value: " << discrete_solution[i] << ". Recursive value: " << "iterative_solution[i]" << std::fixed << std::endl << std::endl;
            correct = false;
        }

    }    
    if(correct) 
        std::cout << "Algorithm completed successfully" << std::endl;
    else 
        std::cout << "Something is wrong!" << std::endl;
}

void
FFT_1D::output_and_test(){
    std::cout << "=================================" << std::endl;
    std::cout << "Printing output:" << std::endl;

    discrete_solve();

    std::cout << "FFT 1D: " << std::endl;
    for (std::size_t i=0; i<N; i++)
    {
        std::cout << "At index " << i << ": " << std::endl;
        std::cout << "Parallel solution: " << parallel_solution[i].real() << " " << parallel_solution[i].imag() << std::endl;
        std::cout << "Discrete solution: " << discrete_solution[i].real() << " " << discrete_solution[i].imag() << std::endl << std::endl;
    }

    test();
}


void
FFT_1D::evaluate_time_and_error(){
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
    std::cout << "Error Evaluation" <<std::endl << "Inverse fft of the parallel solution compared with the initial input " << std::endl;

    cVector inverse_solution = inverse_solve(parallel_solution);
    double max_error = 0;
    for (std::size_t i=0; i<N; i++)
    {
        double error = std::abs(inverse_solution[i] - input[i]);
        if(error > max_error) max_error = error;
    }    
    std::cout << "Max error among all the elements: " << max_error << std::endl;

    //TODO: COMPLETE ERROR COMPUTATION
    //double relative_error = std::abs(parallel_solution. - iterative_solution)/std::abs(iterative_solution);

    //std::cout << "Relative error: " << relative_error << std::endl;
}

void
FFT_1D::save_output_in_file(){
    std::cout << "Writing in output in file output.txt" << std::endl
    std::ofstream file;
    file.open("output.txt");

    for(size_t i=0; i < N; i++){
        file << parallel_solution[i].real() << " " << parallel_solution[i].imag() << std::endl;
    }

    file.close();
}


