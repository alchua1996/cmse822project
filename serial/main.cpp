#include<iostream>
#include<complex>
using std::complex;
using std::string;
typedef complex<double> cplx;
#include <fstream>
using std::ofstream;
#include <chrono>
#include <iomanip>
#include "load_pgm.h"
#include "creatematrix.h"
#include "filters.h"

int main(int argc, char** argv){
    int scales = atoi(argv[1]);
    int orients = atoi(argv[2]);
    string file = argv[3];

    //for rows/cols/count number of filters
    int numrows = 0;
    int numcols = 0;
    
    //stream to load image data
    ofstream outputfile;

    //get image 
    double* img = read_pgm(file, numrows,  numcols);
    
    //start timing from here
    auto start = std::chrono::steady_clock::now();
    
    //allocate a filter bank
    double** filters = new double* [scales*orients+2];
    for(int i = 0; i < scales*orients+2; i++)
    {
        filters[i] = new double[numrows*numcols];
    }
    
    //allocate memory for low and high pass filters in intermediate steps           
    double* highpass = new double[numrows*numcols];
    double* directional = new double[numrows*numcols];
    
    //make low pass at scale J
    filters[0] = create_lowpass(numrows, numcols, scales);

    //make the high pass at scale J
    filters[1] = create_highpass(numrows, numcols, scales);
    
    //make directional filters
    for(int i = 0; i < scales*orients; i++)
    {        
         int j = i/4;
         int q = i%4; 
         filters[i+2] =  create_lowpass(numrows, numcols, j);
         highpass = create_highpass(numrows, numcols, j+1);
         directional = create_directional(numrows, numcols, orients, q);
         filters[i+2] = hadamard(filters[i+2], highpass, numrows, numcols);
         filters[i+2] = hadamard(filters[i+2], directional, numrows, numcols);
    }

    //allocate matrix to store wavelet transform
    double** WT = create_matrix(scales*orients+2, numrows*numcols);
  
    //lowpass part 
    WT[0] =  convolution(filters[0], img, numrows);
    
    //highpass part 
    WT[1] = convolution(filters[1], img, numrows);

    //directional part 
    for(int i =2; i < scales * orients+2; i++)
    {
        WT[i] = convolution(filters[i], img, numrows);

    }

    //take from chrono reference of  https://en.cppreference.com/w/cpp/chrono
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    
    //store this in a file now 
    outputfile.open("WT.txt");
    for(int i = 0; i < scales*orients+2; i++)
    {
        for(int j = 0; j < numrows*numcols; j++)
        {
            outputfile << WT[i][j] << " "; 
        }
        outputfile <<std::endl;
    }
    outputfile.close();
    return 0;
}