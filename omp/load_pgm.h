#include <iostream> // cout, cerr
using std::cerr;
using std::cout;
using std::endl;
#include<string>
using std::string;
#include <fstream> // ifstream
using std::ifstream;
#include <sstream> // stringstream
using std::stringstream;

//NOT gonna lie, a lot of this was from stackexchange. It's just for loading files. Modified for my use.
//https://stackoverflow.com/questions/8126815/how-to-read-in-data-from-a-pgm-file-in-c
double* read_pgm(string file, int& numrows, int& numcols) {
    /*
    Reads in matrix from pgm file
    
    PARAMS: 
    -------
    file: string with filename of pgm file

    OUTPUT:
    -------
    array: pixel values of image saved in a N*M array, where 
    N and M are the number of rows and columns, respectively

    rowNum: saved as pointer in code 

    colNum: saved as pointer in code
    */

    int row = 0;
    int col = 0;
    ifstream infile(file);
    stringstream ss;
    string inputLine = "";

    // First line : version
    getline(infile,inputLine);
    getline(infile,inputLine);

    ss << infile.rdbuf();
    //get rows and cols
    ss >> numcols >> numrows;
    cout << numcols << " columns and " << numrows << " rows" << endl;
    
    double* array =  new double [numcols*numrows];

    //data in file
    for(row = 0; row < numrows; row++)
    {
        for (col = 0; col < numcols; col++)
        {
            ss >> array[row * numcols + col];
        } 
    }
    infile.close();

    //normalize matrix to be between 0 and 1 for stability
    for(row = 0; row < numrows; row++)
    {
        for (col = 0; col < numcols; col++)
        {
            array[row * numcols + col] = static_cast<double>(array[row * numcols + col])/static_cast<double>(255.0);
        } 
    }
    return array;
}