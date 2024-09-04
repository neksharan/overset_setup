#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <cstring>
#include <vector>

const int RECORD_DELIMITER_LENGTH = 4;
const double pi = M_PI;

class Point
{
	public:
		double x, y, z;
		int i, j, k;
};

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

		std::vector<std::vector<int>> acceptor;

        void set_grids(int n,int Npts_x,int Npts_y,const double xmin,const double xmax,const double ymin,const double ymax);
};

class HoleCut 
{
	public:
		double xmin, xmax, ymin, ymax, zmin, zmax;

};


void Read_Grid(Grid*& gl, int& num_blocks, const std::string filename);
void Write_Grid(const double* xs, const double* ys, const double* zs,
                const double* xe, const double* ye, const double* ze,
                const int* Nx_pts, const int* Ny_pts, const int* Nz_pts, 
				const int Nblocks, int base_grid_index, HoleCut& hole_cut, Grid* grid);

void FindFringePoints(Grid*& gl, int num_blocks, int base_grid_index, int fringe_width);
void write_vtk_file(Grid*& gl, const int id, const std::vector<std::vector<int>>& acceptor_pts, std::string base_filename);
void write_points_vtk(std::vector<Point>& points, std::string filename);
// Function to add a vector to the list only if it hasn't been added before
void add_unique(const std::vector<int>& new_point, std::vector<std::vector<int>>& vec);

void Deallocate_grid_memory(Grid*);
