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
	
	Grid* grid = new Grid [Nblocks];

	Write_Grid(xs, ys, zs, 
			   xe, ye, ze, 
			   Nx_pts, Ny_pts, Nz_pts, 
			   Nblocks, grid);	
	Read_Grid("cyl_acoustic_grid.xyz");
	
	return 0;
}
