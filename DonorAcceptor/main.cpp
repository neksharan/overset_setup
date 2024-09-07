#include <Grid.H>
using namespace std;

int main()
{

	int Nblocks = 2;

    int Nx_pts[] = {81, 101};
    int Ny_pts[] = {101, 101};
    int Nz_pts[] = {1, 1};

    double xs[] = {0.5, -10.0};
    double xe[] = {4.0, 10.0};
    double ys[] = {0.0, -10.0};
    double ye[] = {2.0*pi, 10.0};
    double zs[] = {0.0, 0.0};
    double ze[] = {0.0, 0.0};

	int base_grid_index = 1;
	HoleCut hole_cut;
	hole_cut.xmin = -2.5;
	hole_cut.xmax =  2.5;
	hole_cut.ymin = -2.5;
	hole_cut.ymax =  2.5;
	
	Grid* grid_write = new Grid [Nblocks];
	Write_Grid(xs, ys, zs, 
			   xe, ye, ze, 
			   Nx_pts, Ny_pts, Nz_pts, 
			   Nblocks, base_grid_index, hole_cut, grid_write);	

	Grid* grid_read;
	int num_blocks;
	Read_Grid(grid_read, num_blocks,"cyl_acoustic_grid.xyz");

	std::cout  << "num_blocks is " << num_blocks << "\n";

	ComputeCentroids(grid_read, num_blocks);

	// Make alist of fringe points on both grids - ie points which need interpolation 
		
	int fringe_width = 2;
	FindAcceptorPoints(grid_read, num_blocks, base_grid_index, fringe_width);

	//FindDonorPoints(grid_read, num_blocks, base_grid_index);

	return 0;
}
