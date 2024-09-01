#include <DataStruct.H>
#include <cstring>
using namespace std;

double** Create2DMatrix(int Ni, int Nj)
{
    double **the_array = new double* [Ni];
    double *tempxy = new double[Ni*Nj];
 		for ( int i = 0 ; i < Ni; ++i, tempxy += Nj ) {
			the_array[i] = tempxy;
    }

    for(int i(0); i < Ni; ++i)
        for(int j(0); j < Nj; ++j)
		{	the_array[i][j]= 0.;		}

    /*double** the_array = new double* [Ni];
    for(int i(0); i < Ni; ++i)
    {
        the_array[i] = new double[Nj];

        for(int j(0); j < Nj; ++j)
        {
            the_array[i][j] = 0;            
        }
    }*/

    return the_array;
}

int** Create2DMatrix_INT(int Ni, int Nj)
{
    int **the_array = new int* [Ni];
    int *tempxy = new int[Ni*Nj];
 		for ( int i = 0 ; i < Ni; ++i, tempxy += Nj ) {
			the_array[i] = tempxy;
    }

    for(int i(0); i < Ni; ++i)
        for(int j(0); j < Nj; ++j)
		{	the_array[i][j]= 0;		}

    /*int** the_array = new int* [Ni];
    for(int i(0); i < Ni; ++i)
    {
        the_array[i] = new int[Nj];

        for(int j(0); j < Nj; ++j)
        {
            the_array[i][j] = 0;            
        }
    }*/

    return the_array;
}

double*** Create3DMatrix(int Nk, int Ni, int Nj)
{
    double ***the_array = new double**[Nk];
    double **tempxy = new double*[Nk*Ni];
    double *tempxyz = new double[Nk*Ni*Nj];
    for ( int k = 0 ; k < Nk ; ++k, tempxy += Ni ) {
        the_array[k] = tempxy;
		for ( int i = 0 ; i < Ni; ++i, tempxyz += Nj ) {
			the_array[k][i] = tempxyz;
    } }

    for(int k(0); k < Nk; ++k)
        for(int i(0); i < Ni; ++i)
            for(int j(0); j < Nj; ++j)
			{	the_array[k][i][j]= 0.;		}

    return the_array;
}

double**** Create4DMatrix(int Nl, int Ni, int Nj, int Nk)
{
    double ****the_array = new double***[Nl];
    double ***tempxy = new double**[Nl*Ni];
    double **tempxyz = new double*[Nl*Ni*Nj];
    double *tempxyzl = new double[Nl*Ni*Nj*Nk];
    for ( int l = 0 ; l < Nl ; ++l, tempxy += Ni ) {
        the_array[l] = tempxy;
		for ( int i = 0 ; i < Ni; ++i, tempxyz += Nj ) {
			the_array[l][i] = tempxyz;
			for ( int j = 0 ; j < Nj; ++j, tempxyzl += Nk ) {
				the_array[l][i][j] = tempxyzl;
    } } }

    for(int l(0); l < Nl; ++l)
	    for(int i(0); i < Ni; ++i)
	        for(int j(0); j < Nj; ++j)
			    for(int k(0); k < Nk; ++k)
				{	the_array[l][i][j][k]= 0.;		}

    return the_array;
}

int*** Create3DMatrix_INT(int Nk, int Ni, int Nj)
{
    int ***the_array = new int**[Nk];
    int **tempxy = new int*[Nk*Ni];
    int *tempxyz = new int[Nk*Ni*Nj];
    for ( int k = 0 ; k < Nk ; ++k, tempxy += Ni ) {
        the_array[k] = tempxy;
		for ( int i = 0 ; i < Ni; ++i, tempxyz += Nj ) {
			the_array[k][i] = tempxyz;
    } }

    for(int k(0); k < Nk; ++k)
        for(int i(0); i < Ni; ++i)
            for(int j(0); j < Nj; ++j)
			{	the_array[k][i][j]= 0.;		}

    return the_array;
}

void Delete4DMatrix(double****& the_array, int Nl, int Ni, int Nj, int Nk)
{
    delete [] the_array[0][0][0];
    delete [] the_array[0][0];
    delete [] the_array[0];
    delete [] the_array;
}

void Delete3DMatrix(double***& the_array, int Nk, int Ni, int Nj)
{
    delete [] the_array[0][0];
    delete [] the_array[0];
    delete [] the_array;
}

void Delete3DMatrix_INT(int***& the_array, int Nk, int Ni, int Nj)
{
    delete [] the_array[0][0];
    delete [] the_array[0];
    delete [] the_array;
}

void Delete2DMatrix(double**& the_array, int Ni, int Nj)
{
    delete [] the_array[0];
    delete [] the_array;
}

// int and float
static void swap4(void *v)
{
    char    in[4], out[4];
    memcpy(in, v, 4);
    out[0] = in[3];
    out[1] = in[2];
    out[2] = in[1];
    out[3] = in[0];
    memcpy(v, out, 4);
}

// double
static void swap8(void *v)
{
    char    in[8], out[8];
    memcpy(in, v, 8);
    out[0] = in[7];
    out[1] = in[6];
    out[2] = in[5];
    out[3] = in[4];
    out[4] = in[3];
    out[5] = in[2];
    out[6] = in[1];
    out[7] = in[0];
    memcpy(v, out, 8);
}

