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

void Read_Grid(std::string filename);
void Write_Grid(double* xs, double* ys, double* zs,
                double* xe, double* ye, double* ze,
                int* Nx_pts, int* Ny_pts, int* Nz_pts, int Nblocks, Grid* grid);

void Deallocate_grid_memory(Grid*);
