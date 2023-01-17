# FAST FOURIER TRASFORM
## Authors

## **INTRODUCTION**
Here there is a little summarize of our work for the implementation for the **Fast Fourier Trasform**. This Repo includes two folders, one for the 1 dimension fft, and another one for the 2 dimensions fft.
Here below there is a recap of the various files and classes in the folders.

## **FFT_1D**
The folder includes a cpp file containing the class that implements the fft in 1 dimension, an header file for this class and a main file for testing.
In the class called FFT_1D we have different functions and variables to deal with the serial implementation and the parallel one.
Here there's a list of the main ones:

### **Discrete solve**
This fuction takes a vector from the input and implements the discrete fourier trasform on it. It is only use for testing the solution for small input since it's computational cost is high

### **Iterative solve**
Here there's an iterative implemantion of the fast fourier trasform based on the **Cooley-Tuck** algorithm. It is perfomed on a vector saved on the **input** protected variable, and the result it's saved in the **iterative_solution** variable.
This fuction is mainly used to compare the time taken by the serial algorithm, with the time taken by the parallel one.

### **Parallel solve**
This is the parallel implementation of the **Cooley-Tuck** algorithm for the fft. The function uses MPI to scatter the input in samller vectors and distribute them between the processors fo the computer. On every local vector it performs the iterative fft and then gathers the result in the **parallel_solution** variable.
With large-size vectors it performs better then the serial implementation. 

### **Output and test** and **Evaluate time and error**
This two final functions are mainly used to the our algorithm. **Output_and_test** compares the result of the parallel algorithm with the result of the discrete fourier trasform (It shoul be used only with small input due to the computational cost). While the **Evaluate_time_adn_error" calls the iterative solves and compares the time taken by the two algorithms.


### **Compilation**
Finally here there are the instructions to compile the code. In the main file the FFT_1D class is created, then it's loaded a random input and finally it is perfomed the parallel algorithm, followed by the testing function

## **FFT_2D**
In this folder there are the files for the implementation of the class for the FFT in 2 dimensions. In this case we decided to use the **Eigen** library in odrder to deal in an efficient way  with matrixes. Here there's the library documentation: https://eigen.tuxfamily.org/dox/.
the structer of the folder is similar to the one for the 1 dimension fft, since there is an header file for the class, a file containing its implemantation and the main file where the class is created and tested. Most functions uses the function created for the 1D FFT, but are changed in order to work with Eigen matrixes, such as **iterative_solve_wrapped** and **parallel_solve_wrapped**.
Here below there's a little summary of the new functions created:
### **iterative_solve**

### **parallel_solve**

### **Compilation**

## **IMAGE_COMPRESSION**

###




