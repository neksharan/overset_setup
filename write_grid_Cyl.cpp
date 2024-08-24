#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
using namespace std;

const int RECORD_DELIMITER_LENGTH = 4;
const double pi = M_PI;

class Grid
{
	public:
		int Nx, Ny, Nz;
		double*** x;
		double*** y;
		double*** z;
		int*** iblank;
		void Iblanking(int nLev);

		double** rho	;
		double** ux;
		double** uy;
		double** p;
		double** E;
		double*** U;

		void set_grids(int n,int Npts_x,int Npts_y,const double xmin[],const double xmax[],const double ymin[],const double ymax[]);
};

void Deallocate_grid_memory(Grid*);

double** Create2DMatrix(int, int);
int** Create2DMatrix_INT(int, int);
double*** Create3DMatrix(int, int, int);
int*** Create3DMatrix_INT(int, int, int);
double**** Create4DMatrix(int, int, int, int);
void Delete4DMatrix(double****&, int, int, int, int);
void Delete3DMatrix(double***&, int, int, int);
void Delete3DMatrix_INT(int***&, int, int, int);
void Delete2DMatrix(double**&, int, int);

int main()
{
	int Nblocks = 2;
	double dummy, t, r;
	float tt;
	int count = 0;

	Grid* grid = new Grid [Nblocks];

	//int Nx_pts[] = {400, 100};
	//int Ny_pts[] = {300, 100};
	//int Nz_pts[] = {1, 1};
	int Nx_pts[] = {81, 101};
	int Ny_pts[] = {101, 101};
	int Nz_pts[] = {1, 1};
	for (int i=0;i<Nblocks;i++)
	{
		grid[i].Nx = Nx_pts[i];
		grid[i].Ny = Ny_pts[i];
		grid[i].Nz = Nz_pts[i];

		cout << grid[i].Nx << endl;
		cout << grid[i].Ny << endl;
		cout << grid[i].Nz << endl;

		grid[i].x = Create3DMatrix(grid[i].Nx,grid[i].Ny,grid[i].Nz);
		grid[i].y = Create3DMatrix(grid[i].Nx,grid[i].Ny,grid[i].Nz);
		grid[i].z = Create3DMatrix(grid[i].Nx,grid[i].Ny,grid[i].Nz);
		grid[i].iblank = Create3DMatrix_INT(grid[i].Nx,grid[i].Ny,grid[i].Nz);
	}

	//double xs[] = {1.0, -20.0};
	//double xe[] = {18.0, 20.0};
	//double ys[] = {0.0, -20.0};
	//double ye[] = {2.0*pi, 20.0};
	double xs[] = {0.5, -10.0};
	double xe[] = {4.0, 10.0};
	double ys[] = {0.0, -10.0};
	double ye[] = {2.0*pi, 10.0};
	double zs[] = {0.0, 0.0};
	double ze[] = {0.0, 0.0};
	for (int ii=0;ii<Nblocks;ii++)
	{
		if (ii==0)
		{
			double r[grid[ii].Nx];
			for (int i=0;i<grid[ii].Nx;++i)
			{	r[i] = xs[ii] + i * (xe[ii]-xs[ii])/(grid[ii].Nx-1);	}

			double theta[grid[ii].Ny];
			for (int j=0;j<grid[ii].Ny;++j)
			{	theta[j] = ys[ii] + j * (ye[ii]-ys[ii])/(grid[ii].Ny-1);	}

			for (int k=0;k<grid[ii].Nz;++k)
				for (int j=0;j<grid[ii].Ny;++j)
					for (int i=0;i<grid[ii].Nx;++i)
					{				
						grid[ii].x[i][j][k] = r[i]*cos(theta[j]);			
						grid[ii].y[i][j][k] = r[i]*sin(theta[j]);
						grid[ii].z[i][j][k] = 0.0;
						grid[ii].iblank[i][j][k] = 1;
					}
		}
		else
		{
			double dx = (xe[ii]-xs[ii])/(Nx_pts[ii]-1);
			double dy = (ye[ii]-ys[ii])/(Ny_pts[ii]-1);
			double dz = 0.0;
			for (int k=0;k<grid[ii].Nz;++k)
				for (int j=0;j<grid[ii].Ny;++j)
					for (int i=0;i<grid[ii].Nx;++i)
					{
						grid[ii].x[i][j][k] = xs[ii] + dx*i;	
						grid[ii].y[i][j][k] = ys[ii] + dy*j;		
						grid[ii].z[i][j][k] = zs[ii] + dz*k;
						grid[ii].iblank[i][j][k] = 1;
					}
		}
	}

	ofstream outFile;
	char filename[] = "cyl_acoustic_grid.xyz";
	outFile.open (filename, ios::out | ios::binary);

	int size = 4;
	outFile.write((char*) &size, sizeof(int));
	outFile.write((char*) &Nblocks, sizeof(Nblocks));
	outFile.write((char*) &size, sizeof(int));

	size = 4*Nblocks*3;
	outFile.write((char*) &size, sizeof(int));
	for (int i=0;i<Nblocks;i++)
	{			
		outFile.write((char*) &grid[i].Nx, sizeof(grid[i].Nx));
		outFile.write((char*) &grid[i].Ny, sizeof(grid[i].Ny));
		outFile.write((char*) &grid[i].Nz, sizeof(grid[i].Nz));
	}
	outFile.write((char*) &size, sizeof(int));

	for (int ii=0;ii<Nblocks;++ii)
	{
		size = 8*3*grid[ii].Nx*grid[ii].Ny*grid[ii].Nz + 4*grid[ii].Nx*grid[ii].Ny*grid[ii].Nz;
		outFile.write((char*) &size, sizeof(int));
		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
				{
					outFile.write((char*) &grid[ii].x[i][j][k], sizeof(grid[ii].x[i][j][k]));
				}

		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
				{
					outFile.write((char*) &grid[ii].y[i][j][k], sizeof(grid[ii].y[i][j][k]));
				}

		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
				{
					outFile.write((char*) &grid[ii].z[i][j][k], sizeof(grid[ii].z[i][j][k]));
				}

		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
				{
					if (ii>0 && grid[ii].x[i][j][k]>=-2.5 && grid[ii].x[i][j][k]<=2.5 && grid[ii].y[i][j][k]>=-2.5 && grid[ii].y[i][j][k]<=2.5)
					{	grid[ii].iblank[i][j][k] = 0;		}
					outFile.write((char*) &grid[ii].iblank[i][j][k], sizeof(int));
				}

		outFile.write((char*) &size, sizeof(int));
	}

	outFile.close();
	cout << "===============================" << endl;
	printf("Wrote \"%s\" \n", filename);
	cout << "===============================" << endl;
	Deallocate_grid_memory(grid);
	return 0;
}

void Deallocate_grid_memory(Grid* gl)
{
	int Nblocks = 2;
	for (int i=0;i<Nblocks;++i)
	{
		Delete3DMatrix(gl[i].x, gl[i].Nx, gl[i].Ny, gl[i].Nz);
		Delete3DMatrix(gl[i].y, gl[i].Nx, gl[i].Ny, gl[i].Nz);
		Delete3DMatrix(gl[i].z, gl[i].Nx, gl[i].Ny, gl[i].Nz);
		Delete3DMatrix_INT(gl[i].iblank, gl[i].Nx, gl[i].Ny, gl[i].Nz);
	}
}


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


