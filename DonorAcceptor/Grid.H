#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <cstring>
#include <vector>
#include <algorithm>

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
		double*** xcentroid;
		double*** ycentroid;
		double*** zcentroid;
        void Iblanking(int nLev);

        double** rho    ;
        double** ux;
        double** uy;
        double** p;
        double** E;
        double*** U;

		std::vector<std::array<int,3>> acceptors;
		std::vector<std::array<int,3>> donors;

        void set_grids(int n,int Npts_x,int Npts_y,const double xmin,const double xmax,const double ymin,const double ymax);
};

class HoleCut 
{
	public:
		double xmin, xmax, ymin, ymax, zmin, zmax;

};

void Write_Grid(const double* xs, const double* ys, const double* zs,
                const double* xe, const double* ye, const double* ze,
                const int* Nx_pts, const int* Ny_pts, const int* Nz_pts, 
				const int Nblocks, int base_grid_index, HoleCut& hole_cut, Grid* grid);
void Read_Grid(Grid*& gl, int& num_blocks, const std::string filename);
void ComputeCentroids(Grid*& gl, int num_blocks);

void FindAcceptorPoints(Grid*& gl, int num_blocks, int base_grid_index, int fringe_width);
void FindDonorPoints(Grid*& gl, int num_blocks, int base_grid_index);
bool FindDonorsForSinglePoint(Grid*& gl, const int acceptor_grid_id, const int donor_grid_id,
                              const std::array<int,3>& acceptor_pt, std::array<int,3>& idonor);
bool does_cell_enclose_point(const double xa,const double ya, const double za,
                             const Grid& gdon, const std::array<int,3>& cell_to_check);
void write_vtk_file(Grid*& gl, const int id, const std::vector<std::array<int,3>>& acceptor_pts, std::string base_filename);
void write_points_vtk(std::vector<std::array<double,3>>& coords, std::string filename);
// Function to add a vector to the list only if it hasn't been added before
template<typename T1, typename T2>
inline void add_unique(const T1& new_point, T2& vec) {
	 // Check if the vector new_point already exists in vec
    if (std::find(vec.begin(), vec.end(), new_point) == vec.end()) {
        vec.push_back(new_point); // Add the vector if it's not found
    }
}

inline std::array<double,3> ComputeVector(const double x1, const double y1, const double z1,
                                   const double x2, const double y2, const double z2)
{
    double comp0 = x2 - x1;
    double comp1 = y2 - y1;
    double comp2 = z2 - z1;

    double norm = 1.0;//std::sqrt(comp0*comp0 + comp1*comp1 + comp2*comp2);

    std::array<double,3> vec;
    vec[0] = comp0/norm;
    vec[1] = comp1/norm;
    vec[2] = comp2/norm;

    return vec;
}

inline bool is_point_in_box(double x, double y, double z,
                const double xmin, const double xmax, const double ymin,
                const double ymax, const double zmin, const double zmax)
{

    // Shidt the point by a but
    x = x + 1e-10;
    y = y + 1e-10;
    if((x>=xmin and x<=xmax) and
       (y>=ymin and y<=ymax) and
       (z>=zmin and z<=zmax)) {
        return true;
    }

    return false;
}


void Deallocate_grid_memory(Grid*);
