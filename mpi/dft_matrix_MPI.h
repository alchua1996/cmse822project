#include<iostream>
#include<complex>
#include "omp.h"
using std::complex;
typedef complex<double> cplx;

double* transpose(double* matrix, int numRows, int numCols)
{
     
    /*
    Returns transpose of matrix
 
    PARAMS:
    -------
    numRows: number of rows
    numCols: number of columns
    matrix: matrix stored as an numRows*numCols array
    
    OUTPUT:
    -------
    matrix_T: transpose of the matrix
    */
 
    double* matrix_T = new double [numCols*numRows];
    #pragma omp parallel for
    for(int i = 0; i < numRows; i++)
    {     
        for(int j = 0; j <= i; j++)
        {
            if(i != j)
            {
                matrix_T[j * numRows +i] = matrix[i * numCols +j];
                matrix_T[i * numCols +j] = matrix[j * numRows +i];
            }
            else
            {
                matrix_T[i * numCols +j] = matrix[i * numCols +j];
            }
        }
    }
    return matrix_T;
    
}

cplx* transpose(cplx* matrix, int numRows, int numCols)
{
    /*
    Returns transpose of matrix
 
    PARAMS:
    -------
    numRows: number of rows
    numCols: number of columns
    matrix: matrix stored as an numRows*numCols array
    
    OUTPUT:
    -------
    matrix_T: transpose of the matrix
    */
    cplx* matrix_T = new cplx [numCols*numRows];
    #pragma omp parallel for 
    for(int i = 0; i < numRows; i++)
    {   
        for(int j = 0; j <= i; j++)
        {
            if(i != j)
            {
                matrix_T[j * numRows +i] = matrix[i * numCols +j];
                matrix_T[i * numCols +j] = matrix[j * numRows +i];
            }
            else
            {
                matrix_T[i * numCols +j] = matrix[i * numCols +j];
            }
        }
    }
    return matrix_T;
}

double* real(cplx* matrix, int numRows, int numCols)
{
    /*
    Returns real part of matrix, meant to work for compex matrices
 
    PARAMS:
    -------
    numRows: number of rows
    numCols: number of columns
    matrix: matrix stored as an numRows*numCols array
    
    OUTPUT:
    -------
    matrix_real: conjugate of the matrix
    */
 
    double* matrix_real = new double [numCols*numRows];
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < numRows; i++)
    {
        for(int j = 0; j < numCols; j++)
        {   
            matrix_real[i * numCols +j] = real(matrix[i * numCols +j]);
        }
    }
    return matrix_real;
}

double* imag(cplx* matrix, int numRows, int numCols)
{
    /*
    Returns imag part of matrix, meant to work for compex matrices
 
    PARAMS:
    -------
    numRows: number of rows
    numCols: number of columns
    matrix: matrix stored as an numRows*numCols array
    
    OUTPUT:
    -------
    matrix_imag: conjugate of the matrix
    */
 
    double* matrix_imag = new double [numCols*numRows];
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < numRows; i++)
    {
        for(int j = 0; j < numCols; j++)
        {   
            matrix_imag[i * numCols +j] = imag(matrix[i * numCols +j]);
        }
    }
    return matrix_imag;
}


cplx* conj(cplx* matrix, int numRows, int numCols)
{
    /*
    Returns elementwise conjugate of matrix, meant to work for compex matrices
 
    PARAMS:
    -------
    numRows: number of rows
    numCols: number of columns
    matrix: matrix stored as an numRows*numCols array
    
    OUTPUT:
    -------
    matrix_conj: conjugate of the matrix
    */
 
    cplx* matrix_conj = new cplx [numCols*numRows];
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < numRows; i++)
    {
        for(int j = 0; j < numCols; j++)
        {   
            matrix_conj[i * numCols +j] = conj(matrix[i * numCols +j]);
        }
    }
    return matrix_conj;
}

cplx* conj_transpose(cplx* matrix, int numRows, int numCols)
{
    /*
    Returns conjugate transpose of matrix, meant to work for compex matrices
 
    PARAMS:
    -------
    numRows: number of rows
    numCols: number of columns
    matrix: matrix stored as an numRows*numCols array
    
    OUTPUT:
    -------
    matrix_T: conjugate transpose of the matrix
    */
 
    cplx* matrix_T = new cplx [numCols*numRows];
    #pragma omp parallel for 
    for(int i = 0; i < numRows; i++)
    {        
        for(int j = 0; j <= i; j++)
        {   
            if(i != j)
            {
                matrix_T[j * numRows +i] = conj(matrix[i * numCols +j]);
                matrix_T[i * numCols +j] = conj(matrix[j * numRows +i]);
            }
            else
            {
                matrix_T[i * numCols +j] = conj(matrix[i * numCols +j]);
            }
        }
    }
    return matrix_T;
}

double* hadamard(double* a, double* b, int numRows, int numCols)
{
    /*
    Returns hadamard product of two matrices a and b
 
    PARAMS:
    -------
    numRows: number of rows
    numCols: number of columns
    a: matrix stored as an numRows*numCols array
    b: matrix stored as an numRows*numCols array
    
    OUTPUT:
    -------
    result: elementwise product of a and b
    */
    double* result = new double [numRows*numCols];
    #pragma omp parallel for 
    for(int n = 0; n < numRows*numCols; n++)
    {
        result[n] = a[n] * b[n];
    }
    return result;
}

cplx* hadamard(cplx* a, cplx* b, int numRows, int numCols){
    /*
    Returns hadamard product of two matrices a and b.
    This one works for complex doubles instead of normal doubles
 
    PARAMS:
    -------
    numRows: number of rows
    numCols: number of columns
    a: matrix stored as an numRows*numCols array
    b: matrix stored as an numRows*numCols array
    
    OUTPUT:
    -------
    result: elementwise product of a and b
    */

    cplx* result = new cplx [numRows*numCols];
    #pragma omp parallel for 
    for(int n = 0; n < numRows*numCols; n++)
    {
        result[n] = a[n] * b[n];
    }
    return result;
}

//overload, first argument is complex, second is real
cplx* hadamard(cplx* a, double* b, int numRows, int numCols){
    /*
    Returns hadamard product of two matrices a and b.
    This one works for the first argument beign a complex double
    and the second argument being a normal double
 
    PARAMS:
    -------
    numRows: number of rows
    numCols: number of columns
    a: matrix stored as an numRows*numCols array
    b: matrix stored as an numRows*numCols array
    
    OUTPUT:
    -------
    result: elementwise product of a and b
    */
    cplx* result = new cplx [numRows*numCols];
    #pragma omp parallel for
    for(int n = 0; n < numRows*numCols; n++)
    {
        result[n] =  a[n] * b[n];
    }
    return result;
}

//overload, first argument is real, second is complex
cplx* hadamard(double* a, cplx* b, int numRows, int numCols){
    /*
    Returns hadamard product of two matrices a and b.
    This one works for the first argument beign a complex double
    and the second argument being a normal double
 
    PARAMS:
    -------
    numRows: number of rows
    numCols: number of columns
    a: matrix stored as an numRows*numCols array
    b: matrix stored as an numRows*numCols array
    
    OUTPUT:
    -------
    result: elementwise product of a and b
    */    
    cplx* result = new cplx [numRows*numCols];
    #pragma omp parallel for 
    for(int n = 0; n < numRows*numCols; n++)
    {
        result[n] = a[n] * b[n];
    }
    return result;
}

//assume matrices are square and real
double* matmul(double* A, double* B, int numRows)
{
    /*
    Returns matrix multiplication of two matrices two square matrices A and B
 
    PARAMS:
    -------
    numRows: number of rows/cols
    A: matrix stored as an numRows^2 array
    B: matrix stored as an numRows^2 array
    
    OUTPUT:
    -------
    result: A*B stored as numRows^2 array
    */
    double* result = new double [numRows * numRows];
    #pragma omp parallel for collapse(2)
    for(int n = 0; n < numRows; n++)
    {
        for(int m = 0; m < numRows; m++)
        {   
            double sum = 0.0;
            for (int l = 0; l < numRows; l++)
            {
                sum += A[n * numRows + l] * B[l * numRows + m];  
            }
            result[n*numRows+m] = sum;
        }
    }
    return result;   
}

//overload for when both matrices are complex
cplx* matmul(cplx* A, cplx* B, int numRows)
{
    /*
    Returns matrix multiplication of two matrices two square matrices A and B.
    This works when both matrices are complex doubles.
 
    PARAMS:
    -------
    numRows: number of rows/cols
    A: matrix stored as an numRows^2 array
    B: matrix stored as an numRows^2 array
    
    OUTPUT:
    -------
    result: A*B stored as numRows^2 array
    */
    cplx* result = new cplx [numRows * numRows];
    #pragma omp parallel for collapse(2)
    for(int n = 0; n < numRows; n++)
    {
        for(int m = 0; m < numRows; m++)
        {   
            cplx mysum (0.0,0.0);
            for (int l = 0; l < numRows; l++)
            {
                mysum += A[n * numRows + l] * B[l * numRows + m];      
            }
            result[n*numRows+m] = mysum;
            
        }
    }
    return result;      
}

//overload when only first matrix is complex
cplx* matmul(cplx* A, double* B, int numRows)
{
    /*
    Returns matrix multiplication of two matrices two square matrices A and B.
    This works when both matrices are complex doubles.
 
    PARAMS:
    -------
    numRows: number of rows/cols
    A: matrix stored as an numRows^2 array
    B: matrix stored as an numRows^2 array
    
    OUTPUT:
    -------
    result: A*B stored as numRows^2 array
    */
    cplx* result = new cplx [numRows * numRows];
    #pragma omp parallel for collapse(2)
    for(int n = 0; n < numRows; n++)
    {
        for(int m = 0; m < numRows; m++)
        {   
            cplx mysum (0.0,0.0);
            for (int l = 0; l < numRows; l++)
            {
                mysum += A[n * numRows + l] * B[l * numRows + m];      
            }
            result[n*numRows+m] = mysum;
        }
    }
    return result;        
}

//overload when only second matrix is complex
cplx* matmul(double* A, cplx* B, int numRows)
{
    /*
    Returns matrix multiplication of two matrices two square matrices A and B.
    This works when both matrices are complex doubles.
 
    PARAMS:
    -------
    numRows: number of rows/cols
    A: matrix stored as an numRows^2 array
    B: matrix stored as an numRows^2 array
    
    OUTPUT:
    -------
    result: A*B stored as numRows^2 array
    */
    cplx* result = new cplx [numRows * numRows];
    #pragma omp parallel for collapse(2)
    for(int n = 0; n < numRows; n++)
    {
        for(int m = 0; m < numRows; m++)
        {   
            cplx mysum (0.0,0.0);
            for (int l = 0; l < numRows; l++)
            {
                mysum += A[n * numRows + l] * B[l * numRows + m];  
            }
            result[n*numRows+m] = mysum;
        }
    }
    return result;      
}

//make dft matrix
cplx* DFT_matrix (int numRows)
{
    /*
    Returns DFT matrix for length N signal
 
    PARAMS:
    -------
    numRows: number of rows
    
    OUTPUT:
    -------
    DFT_mat: DFT matrix for signal of length N
    */
 
    //declare matrix
    cplx* DFT_mat;
    //declare constant for i
    const cplx i (0, 1);
    //declare constant for one
    const cplx one (1,0);
    //get the specified root
    const cplx w = std::exp(2 * M_PI * i  * static_cast<double>(1.0)/static_cast<double>(numRows));
 
    //allocate the memory for DFT
    DFT_mat = new cplx [numRows*numRows];
    //put stuff into dft matrix
    #pragma omp parallel for collapse(2)
    for(int n = 0; n < numRows; n++)
    {
        for(int m = 0; m < numRows; m++)
        {
            if(n == 0 or m == 0)
            {
                DFT_mat[n*numRows + m] = one/std::sqrt(numRows);
            }
            else
            {
                DFT_mat[n*numRows + m] = std::pow(w, n*m)/std::sqrt(numRows);
            }
        }    
    }
    return DFT_mat; 
}

//calculates DFT
cplx* DFT(double* signal, int numRows)
{
    /*
    If my code wasn't trash, I would make a function where I don't need to
    reallocate memory each time I make a matrix multiplication, but I'm honestly
    way too lazy at this point. There is probably a way of overloading.
 
    This calculate the DFT using a DFT matrix.
 
    INPUTS:
    -------
    signal: input signal to take DFT of
    numRows: number of rows, could be determined from the signal, but I don't
    really care because I'm running out of time to write this.

 
    OUTPUTS:
    --------
    signal_DFT: complex double* with numRows*numRows entries containing the DFT
    of the signal.
    */
  
    
    //declare array for DFT matrix
    cplx* DFT_mat = DFT_matrix(numRows);

    //do matrix multiplication to calculate DFT on cols first
    cplx* signal_DFT = matmul(DFT_mat, signal, numRows);
 
    //calculate DFT on columns now
    signal_DFT = matmul(signal_DFT, DFT_mat, numRows);
    signal_DFT = conj(signal_DFT, numRows, numRows);
 
    delete [] DFT_mat;
    return signal_DFT;
}

//calculate IDFT
double* IDFT(cplx* signal, int numRows)
{
    /*
    If my code wasn't trash, I would make a function where I don't need to
    reallocate memory each time I make a matrix multiplication, but I'm honestly
    way too lazy at this point. There is probably a way of overloading.
 
    This calculate the IDFT using a IDFT matrix.
 
    INPUTS:
    -------
    signal: input signal to take IDFT of
    numRows: number of rows, could be determined from the signal, but I don't
    really care because I'm running out of time to write this.

 
    OUTPUTS:
    --------
    signal_IDFT: double* with numRows*numRows entries containing the DFT
    of the signal.
    */
  
    
    //declare array for IDFT matrix
    //make dft matrix
    cplx* IDFT_mat = DFT_matrix(numRows);
 
    //calculate conjugate transpose since the operator is unitary
    IDFT_mat = conj_transpose(IDFT_mat, numRows, numRows);
    
    //undo complex conjugate
    cplx*signal_conj = conj(signal, numRows, numRows);
     
    //undo right matrix multiplcation
    cplx* signal_IDFT = matmul(signal_conj, IDFT_mat, numRows);
 
    //undo left matrix multiplication
    signal_IDFT = matmul(IDFT_mat, signal_IDFT, numRows);

    //clear the memory
    delete [] IDFT_mat;
    delete [] signal_conj;
 
    //take the real part because I want a real signal
    return real(signal_IDFT, numRows, numRows);
}

//calculate shift of zero frequency
void shift(cplx* signal, int numRows)
{
    /*
 
    Shifts image so that zero frequency is at center, not edge. Assume
    signal has an even number of elements in each row and column.
 
    INPUTS:
    -------
    filter: signal to be shifted
    numRows: number of rows, could be determined from the signal

 
    OUTPUTS:
    --------
    result after shifting is in place, so no need to return
    */
 
    //get the middle row index
    int middle = numRows/2;
    #pragma omp parallel for collapse(2)
    for(int i = 0;i < middle; i++)
    {
        for(int j = 0; j < middle; j++)
        {
            //swap elements in top left and bottom right
            cplx temp1 = signal[i*numRows + j];
            signal[i*numRows + j] = signal[(i+middle) *numRows + j+middle];
            signal[(i+middle) *numRows + j+middle] = temp1;
         
            //swap elements in top right and bottom left
            cplx temp2 = signal[i *numRows + j+middle];
            signal[i *numRows + j+middle] = signal[(i+middle)*numRows + j];
            signal[(i+middle)*numRows + j] = temp2;           
        }
    }
}

//calculate convolution in frequency now
double* convolution(double* filter, double* signal, int numRows)
{
    /*
 
    Calculates convolution in frequency. 
 
    INPUTS:
    -------
    filter: Heeger Bergen wavelets given in frequency
    signal: input signal to be convolved
    numRows: number of rows, could be determined from the signal

 
    OUTPUTS:
    --------
    result: result from convolution
 
    Method:
    -------
    Take filter in frequency and call it f. Take signal in space and call it g.
    Find the DFT of g, and then fftshift. Calculuate the hadamard product f*g. 
    Then take the ifft and return the result.
    */
 
    //calculate dft of signal
    cplx* signal_DFT = DFT(signal, numRows);
    
    //fftshift now
    shift(signal_DFT, numRows);
 
    //hadamard product
    signal_DFT = hadamard(signal_DFT, filter, numRows, numRows);
   
    //desired result after idft
    double* result = IDFT(signal_DFT, numRows);
   
    //free the memory
    delete [] signal_DFT;
 
    return result;
}