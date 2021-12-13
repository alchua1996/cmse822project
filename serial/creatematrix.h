//I used stackoverflow for reference on how to do this.
//https://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new
double** create_matrix(int row, int col)
{
    /*
    Allocates memory for a 2D matrix 
 
    PARAMS:
    -------
    row: number of rows for matrix
    columns: number of columns for matrix
    
    OUTPUT:
    -------
    matrix: 2D matrix with rows by columns memory allocated
    */
 
    double** matrix = new double* [row]; 
    for (int i = 0; i < row; ++i)
    {
        matrix[i] = new double[col](); 
    }
    return matrix;
}

void destroy_matrix(double** &mat, int rows)
{
    /*
    Frees memory from matrix
 
    PARAMS:
    -------
    mat: 2D matrix
    rows: number of rows in matrix
    
    OUTPUT:
    -------
    matrix is deallocated, no outputs
    */
 
    if (mat)
    {
        for (int i = 0; i < rows; ++i)
        {
            delete[] mat[i]; //delete each row..
        }

        delete[] mat;  //delete the rows..
        mat = nullptr;
    }
}