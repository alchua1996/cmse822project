#include<iostream>
#include<complex>
using std::complex;
using std::string;
typedef complex<double> cplx;
#include <fstream>
using std::ofstream;

#include <iomanip>
#include "load_pgm.h"
#include "creatematrix.h"
#include "filters_MPI.h"
#include "omp.h"
#include <mpi.h>

int main(int argc, char** argv){
    int scales = atoi(argv[1]);
    int orients = atoi(argv[2]);
    int tcount = atoi(argv[3]);
    int numrows = atoi(argv[4]);
    int numcols = atoi(argv[5]);
    string file = argv[6];
    
    //stream to load image data
    ofstream outputfile;

    //MPI calls
    int comm_sz;
	  int myrank;
    MPI_Status status;
	
	  //mpi calls
	  MPI_Init(&argc, &argv);
	  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    //start timing from here
    double start = 0;
    double end = 0;
    
    //declare image pointer
    double* img;
    
    //get image 
    img = read_pgm(file, numrows,  numcols);
    
    //allocate memory for low and high pass filters in intermediate steps           
    double* highpass;
    double* directional;
    
    //filter used to store info during the scatter
    double* filter_scatter;
    //filter used to store wavelet transform during scattoer
    double* WT_scatter;
    
    //allocate storage for wavelet transform
    double* WT;
     
    if(myrank == 0) //make low pass filter and allocate memory in master rank
    {  
        //allocate matrix to store wavelet transform
        WT =  new double[(scales*orients + 2)*numrows*numcols];
        
        //set thread number 
        omp_set_num_threads(tcount);
 
        //start timing
        start = MPI_Wtime();    
    }

     
    //allocate memory
    highpass = new double[numrows*numcols];
    directional = new double[numrows*numcols];
    filter_scatter = new double[numrows*numcols];
    WT_scatter = new double[numrows*numcols];
    
 
     
    if(myrank == 0)
    {
        create_lowpass(filter_scatter, numrows, numcols, scales);
        WT_scatter = convolution(filter_scatter, img, numrows);
    }
    else if(myrank == 1)
    {
        create_highpass(filter_scatter, numrows, numcols, scales);
        WT_scatter = convolution(filter_scatter, img, numrows);
    }
    else if(myrank > 1 && myrank < scales*orients + 2)
    {
        int j = (myrank-2)/4;
        int q = (myrank-2)%4; 
        create_lowpass(filter_scatter, numrows, numcols, j);
        create_highpass(highpass, numrows, numcols, j+1);
        create_directional(directional, numrows, numcols, orients, q);
        filter_scatter = hadamard(filter_scatter, highpass, numrows, numcols);
        filter_scatter = hadamard(filter_scatter, directional, numrows, numcols);
        WT_scatter = convolution(filter_scatter, img, numrows);  
    }
    
    //gather into array
    MPI_Gather(WT_scatter, numrows*numcols, MPI_DOUBLE, WT, numrows*numcols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    //all the processes finish and gather to rank 0 process for other parts
    if(myrank == 0)
    {
        //finish timing
        end = MPI_Wtime();
        auto elapsed_seconds = end-start;
        std::cout << "elapsed time: " << elapsed_seconds << "s\n";
        
        //store results in a file now 
        outputfile.open("WT_test.txt");
        for(int i = 0; i < scales*orients+2; i++)
        {
            for(int j = 0; j < numrows*numcols; j++)
            {
                outputfile << WT[i * numrows*numcols + j] << " "; 
            }
            outputfile <<std::endl;
        }
        outputfile.close();
    
        //clear the memory
        delete [] filter_scatter;
        delete [] WT_scatter;
        delete [] highpass;
        delete [] directional;
    }
    MPI_Finalize();
    return 0;
}