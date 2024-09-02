#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <cstring>

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

        double** rho    ;
        double** ux;
        double** uy;
        double** p;
        double** E;
        double*** U;

        void set_grids(int n,int Npts_x,int Npts_y,const double xmin,const double xmax,const double ymin,const double ymax);
};

class HoleCut 
{
	public:
		double xmin, xmax, ymin, ymax, zmin, zmax;

};

void Read_Grid(Grid* gl, int num_blocks, const std::string filename);
void Write_Grid(const double* xs, const double* ys, const double* zs,
                const double* xe, const double* ye, const double* ze,
                const int* Nx_pts, const int* Ny_pts, const int* Nz_pts, 
				const int Nblocks, int base_grid_index, HoleCut& hole_cut, Grid* grid);

void Deallocate_grid_memory(Grid*);
