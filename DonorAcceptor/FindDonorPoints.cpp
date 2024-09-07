#include <Grid.H>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;


void FindDonorsSinglePoint(Grid*& gl, int acceptor_grid_id, int donor_grid_id, std::array<int,3>& acceptor_pt)
{

	int i = acceptor_pt[0];	
	int j = acceptor_pt[1];	
	int k = acceptor_pt[2];

	double x = gl[acceptor_grid_id].x[i][j][k]; 
	double y = gl[acceptor_grid_id].y[i][j][k]; 
	double z = gl[acceptor_grid_id].z[i][j][k];

	

}



void FindDonorPoints(Grid*& gl, int num_blocks, int base_grid_index)
{

	// Consider the acceptor points on the base grid. For each of these points, find the cell 
	// in the other grid which encloses the acceptor point

	int id = base_grid_index;	

	std::cout << "Size of acceptors is " << gl[id].acceptors.size() << "\n";
	for(int pt=0;pt<gl[id].acceptors.size();pt++){
		std::array<int,3> acceptor_pt = gl[id].acceptors[pt];
		//FindDonorPointsForSinglePoint(gl,1,0,acceptor_pt);				
	}

}


