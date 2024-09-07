#include <Grid.H>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

void write_vtk_file(Grid*& gl, const int id, const std::vector<std::vector<int>>& acceptor_pts, std::string base_filename) {

	std::string filename = base_filename + "_" + std::to_string(id) + ".vtk";
	std::vector<Point> points;
	points.resize(acceptor_pts.size());
	for(int pt=0;pt<acceptor_pts.size();pt++){
		int i = acceptor_pts[pt][0];
		int j = acceptor_pts[pt][1];
		int k = acceptor_pts[pt][2];
		points[pt].x = gl[id].x[i][j][k];
		points[pt].y = gl[id].y[i][j][k];
		points[pt].z = 1e-12;//gl[id].z[i][j][k];
	}	

	write_points_vtk(points, filename.c_str());
}

void FindAcceptorPoints(Grid*& gl, int num_blocks, int base_grid_index, int fringe_width)
{

	int id = base_grid_index;	
	int fw = fringe_width;
	std::cout << "Fringe width is " << fw << "\n";

	// Fund the internal boundary

	int imin = -1;
	int imax = -1;
	int jmin = -1;
	int jmax = -1;
	int kmin = 0;
	int kmax = 0;

	for (int i=1;i<gl[id].Nx-1;++i) {
		for (int j=1;j<gl[id].Ny-1;++j) {
			for (int k=0;k<gl[id].Nz;++k) {	
				if(gl[id].iblank[i][j][k] == 1 and gl[id].iblank[i+1][j][k] == 0) {
					imin = i;		
				}
				if(gl[id].iblank[i][j][k] == 1 and gl[id].iblank[i-1][j][k] == 0) {
					imax = i;
				}
				if(gl[id].iblank[i][j][k] == 1 and gl[id].iblank[i][j+1][k] == 0) {
					jmin = j;		
				}
				if(gl[id].iblank[i][j][k] == 1 and gl[id].iblank[i][j-1][k] == 0) {
					jmax = j;
				}
			}
		}
	}

	std::cout << imin  << " " << imax << " " << jmin << " " << jmax << "\n";

	std::vector<std::vector<int>>& acceptor_pts_0 = gl[0].acceptors;
	std::vector<std::vector<int>>& acceptor_pts_1 = gl[1].acceptors;

	std::vector<int> new_point;
		
	for (int f=0; f<fw; f++) {
		for (int i=imin-f;i<=imax+f;i++) {
			new_point = {i,jmin-f,0};
			add_unique(new_point, acceptor_pts_1); 
			new_point = {i,jmax+f,0};
			add_unique(new_point, acceptor_pts_1);
		}
		for (int j=jmin-f;j<=jmax+f;j++) {
			new_point = {imin-f,j,0};
			add_unique(new_point, acceptor_pts_1);
			new_point = {imax+f,j,0};
			add_unique(new_point, acceptor_pts_1);
		}
	}

	write_vtk_file(gl, id, acceptor_pts_1, "acceptor_pts");

	// On the fine grid, for now just take the last 2 layers of grid points in the first dimensioan
	// This is some hard coding, but alright for now.

	id = 0;

	for (int i=gl[id].Nx-fw;i<gl[id].Nx;++i) {
		for (int j=1;j<gl[id].Ny;++j) {
			for (int k=0;k<gl[id].Nz;++k) {	
				new_point = {i,j,k};
				add_unique(new_point, acceptor_pts_0);
			}
		}
	}
		
	write_vtk_file(gl, id, acceptor_pts_0, "acceptor_pts");	
}


