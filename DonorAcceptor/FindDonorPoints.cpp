#include <Grid.H>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

bool does_cell_enclose_point(const double xa,const double ya, const double za, 
						     const Grid& gdon, const std::array<int,3>& cell_to_check)
{

	// For a point to be inside a cell, the dot product of the vector joining point 
	// to the cell face centers, and the vector joining the cell centroids to the cell
	// face centers should be positive, for all the faces (taking advantage of the 
	// fact that all cells are convex)

	std::array<double,3> norm_check, face_norm;
	std::array<double,3> centroid;
	std::array<double,4> x, y, xcen, ycen;
	
	int i= cell_to_check[0];
	int j= cell_to_check[1];
	int k= cell_to_check[2];

	x[0] = gdon.x[i][j][k]; x[1] = gdon.x[i+1][j][k]; x[2] = gdon.x[i+1][j+1][k]; x[3] = gdon.x[i][j+1][k];
	y[0] = gdon.y[i][j][k]; y[1] = gdon.y[i+1][j][k]; y[2] = gdon.y[i+1][j+1][k]; y[3] = gdon.y[i][j+1][k];

	xcen[0] = 0.5*(x[0] + x[1]); xcen[1] = 0.5*(x[1] + x[2]); xcen[2] = 0.5*(x[2] + x[3]); xcen[3] = 0.5*(x[3] + x[0]);
	ycen[0] = 0.5*(y[0] + y[1]); ycen[1] = 0.5*(y[1] + y[2]); ycen[2] = 0.5*(y[2] + y[3]); ycen[3] = 0.5*(y[3] + y[0]);

	centroid = {gdon.xcentroid[i][j][k],
			    gdon.ycentroid[i][j][k],
				gdon.zcentroid[i][j][k]};


	//printf("Point checking. %0.15f, %0.15f, %0.15f\n",xa,ya,za);
	
	for(int f=0;f<4;f++){
		face_norm  = ComputeVector(centroid[0],centroid[1],0.0,xcen[f],ycen[f],0.0);
		norm_check = ComputeVector(xa, ya, 0.0, xcen[f], ycen[f], 0.0); 	
	
		// Find dot product
		double dot_prod = 0.0;
		for(int dir=0;dir<3;dir++){
			dot_prod = dot_prod + face_norm[dir]*norm_check[dir];
		}
		/*std::cout << "Dot prod with " << f << " " << dot_prod << "\n";
		std::cout << "Face norm " << face_norm[0] << " " << face_norm[1] << "\n";
		std::cout << "norm check " << norm_check[0] << " " << norm_check[1] << "\n";
		if(f<=2){
			std::cout << "coords " << x[f] << " " << x[f+1] << " " << y[f] << " " << y[f+1] << "\n";
		}
		else if(f==3) {
			std::cout << "coords " << x[f] << " " << x[0] << " " << y[f] << " " << y[0] << "\n";
		}
		std::cout << "Face center " << xcen[f] << " " << ycen[f] << "\n";
		std::cout << "Centroids " << centroid[0] << " " << centroid[1] << "\n";*/
		if (dot_prod < 0.0){
			return false;
		}	
	}
	//std::cout << "Cell found "<< "\n";
	return true;		
}

bool FindDonorsForSinglePoint(Grid*& gl, const int acceptor_grid_id, const int donor_grid_id, 
							  const std::array<int,3>& acceptor_pt, std::array<int,3>& idonor)
{

	// Loop over all aceptor points of the grid 

	const Grid &gacc = gl[acceptor_grid_id];
	int ia = acceptor_pt[0];	
	int ja = acceptor_pt[1];	
	int ka = acceptor_pt[2];

	double xa = gacc.x[ia][ja][ka] + 1e-6; 
	double ya = gacc.y[ia][ja][ka] + 1e-6; 
	double za = gacc.z[ia][ja][ka];

	std::vector<std::array<double,3>> point;
	point.resize(1,{xa,ya,za});

	write_points_vtk(point,"point_check.vtk");	
	
	// Find the cell in the other grid which encloses the acceptor point

	double width = 0.5;

	double xmin, xmax, ymin, ymax, zmin, zmax;

	xmin = xa-width/2.0;
	xmax = xa+width/2.0;
	ymin = ya-width/2.0;
	ymax = ya+width/2.0;
	zmin = 0.0;
	zmax = 0.0;
	
	const Grid &gdon = gl[donor_grid_id];
	std::vector<std::array<int,3>> cells_to_check;

	for(int i=0;i<gdon.Nx-1;i++) {
		for(int j=0;j<gdon.Ny-1;j++) {
			for(int k=0;k<gdon.Nz;k++) {
				// Find if the centroid is within the box
				std::array<double,3> centroid = {gdon.xcentroid[i][j][k],
												  gdon.ycentroid[i][j][k],
												  gdon.zcentroid[i][j][k]};
				if(is_point_in_box(centroid[0],centroid[1],centroid[2],
								   xmin,xmax,ymin,ymax,zmin,zmax)) {
					cells_to_check.push_back({i,j,k});
				}
			}
		}
	}	
				
	std::vector<std::array<double,3>> check_donor_cells;
	for(int pt=0;pt<cells_to_check.size();pt++){
		int i = cells_to_check[pt][0];
		int j = cells_to_check[pt][1];
		int k = cells_to_check[pt][2];
		check_donor_cells.push_back({gdon.xcentroid[i][j][k],
									 gdon.ycentroid[i][j][k],	
									 gdon.zcentroid[i][j][k]});
	}
	
	write_points_vtk(check_donor_cells,"donor_cells.vtk");	
	
	// TODO: Do bounding box check

	
	// Do final enclosing check
	check_donor_cells.clear();

	for(int pt=0;pt<cells_to_check.size();pt++){
	//for(int pt=21;pt<22;pt++){
		std::array<int,3>& cell = cells_to_check[pt];
		if(does_cell_enclose_point(xa,ya,za,gdon,cell)) {
			idonor = cell;
			int i = cell[0];
    	    int j = cell[1];
        	int k = cell[2];
			check_donor_cells.push_back({gdon.xcentroid[i][j][k],
									 gdon.ycentroid[i][j][k],	
									 gdon.zcentroid[i][j][k]});
				
	 		write_points_vtk(check_donor_cells,"final_donor_cell.vtk");
			return true;
		}
	}

	return false;
}

void FindDonorPoints(Grid*& gl, int num_blocks, int base_grid_index)
{

	// Consider the acceptor points on the base grid. For each of these points, find the cell 
	// in the other grid which encloses the acceptor point

	int donor_grid_id;
	for(int id=0;id<num_blocks;id++){
		std::cout << "Size of acceptors is " << gl[id].acceptors.size() << "\n";
		std::array<int,3> idonor;
		if(id==0) {
			donor_grid_id = 1;
		}
		else if(id==1) {
			donor_grid_id = 0;
		}
		for(int pt=0;pt<gl[id].acceptors.size();pt++){
		//for(int pt=26;pt<27;pt++){
			std::array<int,3> acceptor_pt = gl[id].acceptors[pt];
			bool has_donors = FindDonorsForSinglePoint(gl,id,donor_grid_id,acceptor_pt,idonor);
			if(has_donors){
				gl[id].donors.push_back(idonor);
			}
			else if(!has_donors){
				std::cout << "Found orphan point " << pt << " " << acceptor_pt[0] << " " << acceptor_pt[1] << " " << acceptor_pt[2] << ". Exiting..." << "\n";
				exit(0);
			}				
		}
		std::vector<std::array<double,3>> donor_pts;
		for(int pt=0;pt<gl[id].donors.size();pt++){
			std::array<int,3> idonor = gl[id].donors[pt];
			int i = idonor[0];
			int j = idonor[1];
			int k = idonor[2];
			donor_pts.push_back({gl[donor_grid_id].xcentroid[i][j][k],
								gl[donor_grid_id].ycentroid[i][j][k],
								0.0});
		}
		std::string basefilename = "donor_pts";
		std::string filename = basefilename + "_" + std::to_string(id) + ".vtk";
		write_points_vtk(donor_pts,filename); 
	}
}


