#include <Grid.H>
#include <DataStruct.H>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

void ComputeCentroids(Grid*& gl, int num_blocks)
{

	for(int id=0;id<num_blocks;id++) {
		gl[id].xcentroid = Create3DMatrix(gl[id].Nx-1,gl[id].Ny-1,gl[id].Nz);	
		gl[id].ycentroid = Create3DMatrix(gl[id].Nx-1,gl[id].Ny-1,gl[id].Nz);	
		gl[id].zcentroid = Create3DMatrix(gl[id].Nx-1,gl[id].Ny-1,gl[id].Nz);	
		std::vector<std::array<double,3>> centroids;
		for (int i=0;i<gl[id].Nx-1;++i) {
			for (int j=0;j<gl[id].Ny-1;++j) {
				for (int k=0;k<gl[id].Nz;++k) {	
					std::array<double,3> centroid;
					centroid[0] = 0.25*(gl[id].x[i][j][k] + gl[id].x[i+1][j][k] +
							  		    gl[id].x[i+1][j+1][k] + gl[id].x[i][j+1][k]);	

					centroid[1] = 0.25*(gl[id].y[i][j][k] + gl[id].y[i+1][j][k] +
							            gl[id].y[i+1][j+1][k] + gl[id].y[i][j+1][k]);	

					centroid[2] = 0.0;
					gl[id].xcentroid[i][j][k] = centroid[0];	
					gl[id].ycentroid[i][j][k] = centroid[1];	
					gl[id].zcentroid[i][j][k] = centroid[2];
					centroids.push_back({centroid[0],centroid[1],1e-12});
				}
			}
		}
		std::string basefilename = "grid_centroids";
		std::string filename = basefilename + "_" + std::to_string(id) + ".vtk"; 
		write_points_vtk(centroids,filename);	
	}
}


