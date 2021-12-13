#define _USE_MATH_DEFINES
#include<iostream>
#include <cmath>
#include <iomanip>
#include "matrix_dft_ops.h"
using std::pow;
using std::abs;
using std::sqrt;
using std::cos;
using std::log2;
using std::atan2;

//why is there no factorial built into C++...
int factorial(int n)
{   
    if(n == 0 or n < 0)
    {
        return 1;
    }
    else
    {
        int output = 1;
        for(int i = 2; i <= n; i ++)
        {
            output *= i;    
        }
        return output; 
    }
}

double* create_lowpass(int rows, int cols, int scale)
{   
    /*
    Creates low pass filter specified in Heeger Bergen paper
    
    PARAMS:
    -------
    rows: number of rows for filter
    cols: number of columns for filter
    scale: scale/zoom of filter
 
    OUTPUT:
    -------
    filter: low-pass filter at given scale with specified rows/cols
    */
 
    //allocate memory
    double* filter = new double [rows*cols];

    //generate filter using double loop 
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++){
            double r = std::sqrt(std::pow(2.0 * M_PI * (i - rows/2)/rows,2) + std::pow(2.0 * M_PI * (j - cols/2)/cols,2));
            double z = 1.0 * pow(2.0,1.0 * scale) * r;
            
            if(z <= M_PI/2)
            {
                filter[i * rows + j] = 1.0;
            }
            else if((M_PI/2 < z)  && (z <= M_PI))
            {
                filter[i * rows + j] = cos(M_PI/2 * log2(2.0 * pow(2,scale) * r/ M_PI) );
            }
            else
            {
                filter[i * rows + j] = 0.0;
            }   
        }
    }
    return filter;    
}


double* create_highpass(int rows, int cols, int scale)
{   
    /*
    Creates high pass filter specified in Heeger Bergen paper
    
    PARAMS:
    -------
    rows: number of rows for filter
    cols: number of columns for filter
    scale: scale/zoom of filter
 
    OUTPUT:
    -------
    filter: high-pass filter at given scale with specified rows/cols
    */
 
    //allocate memory
    double* filter = new double [rows*cols];

    //generate filter using double loop 
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++){
            double r = std::sqrt(std::pow(2.0 * M_PI * (i - rows/2)/rows,2) + std::pow(2.0 * M_PI * (j - cols/2)/cols,2));
            double z = 1.0 * pow(2.0,1.0 * scale) * r;
            
            if(z <= M_PI/2)
            {
                filter[i * rows + j] = 0.0;
            }
            else if((M_PI/2 < z)  && (z <= M_PI))
            {
                filter[i * rows + j] = cos(M_PI/2 * log2(pow(2,scale) * r/ M_PI) );
            }
            else
            {
                filter[i * rows + j] = 1.0;
            } 
        }
    }
    return filter;    
}

double* create_directional(int rows, int cols, int Q, int q)
{
    /*
    Creates high pass filter specified in Heeger Bergen paper
    
    PARAMS:
    -------
    rows: number of rows for filter
    cols: number of columns for filter
    scale: scale/zoom of filter
 
    OUTPUT:
    -------
    filter: high-pass filter at given scale with specified rows/cols
    */
 
    double* filter = new double [rows*cols];
    const double ALPHA = pow(2, Q-1) * factorial(Q-1)/ sqrt(Q*factorial(2*Q-2));
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++){
            //calculate angle
            double x = 2.0 * M_PI * (i - rows/2)/rows;
            double y = 2.0 * M_PI * (j - cols/2)/cols;
            double theta = atan2(y, x);
            
            //create values for indicator function
            double r1 = theta - M_PI * static_cast<double>(q)/(1.0*Q);
            double r2 = theta - M_PI * static_cast<double>(q-Q)/(Q *1.0);

            //ensure values are between -pi and pi
            if(r1 < -M_PI)
            {
                r1 += 2 * M_PI;
            }
            if(r2 > M_PI){
                r2 -= 2 * M_PI;
            }

            //create the filters
            double z1 = ALPHA *pow(cos(r1), Q-1);
            double z2 = ALPHA * pow(cos(r2), Q-1);
            if(abs(r1) <= M_PI/2.0 && abs(r2) <= M_PI/2.0)
            {
                filter[i * cols + j] = z1 + z2;
            }
            else if(abs(r1) <= M_PI/2.0 && abs(r2) > M_PI/2.0)
            {
                filter[i * cols + j] = z1;
            }
            else if(abs(r1) > M_PI/2.0 && abs(r2) <= M_PI/2.0)
            {
                filter[i * cols + j] = z2;
            }
            else
            {
                filter[i * cols + j] = 0.0;
            }   
        }
    }  
    return filter;  
}

double** createPsi(int numRows, int numCols, int J, int Q){
    
    //allocate array for filter bank
    //each row has a filter in it that is numRows*numCols
    //there are J*Q+2 rows in total

    double** filterBank = new double* [J*Q+2];
    for(int i = 0; i < J*Q+2; i++)
    {
        filterBank[i] = new double[numRows*numCols];
    }
    
    //allocate memory for low and high pass filters in intermediate steps           
    double* highpass = new double[numRows*numCols];
    double* directional = new double[numRows*numCols];
    
    //make low pass at scale J
    filterBank[0] = create_lowpass(numRows, numCols, J);

    //make the high pass at scale J
    filterBank[1] = create_highpass(numRows, numCols, J);
    
    //make directional filters
    int count = 2;
    for(int j = 0; j < J; j++)
    {
        for(int q = 0; q < Q; q++)
        {
            
            filterBank[count] =  create_lowpass(numRows, numCols, j);
            double* highpass = create_highpass(numRows, numCols, j+1);
            double* directional = create_directional(numRows, numCols, Q, q);
            filterBank[count] = hadamard(filterBank[count], highpass, numRows, numCols);
            filterBank[count] = hadamard(filterBank[count], directional, numRows, numCols);
            count += 1;
        }
    }
    
    //clear the memory
    delete [] highpass;
    delete [] directional;
    
    //output filterbank
    return filterBank;
}