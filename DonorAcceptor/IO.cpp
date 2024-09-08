#include <Grid.H>
#include <vector>
#include <iostream>
using namespace std;

void write_vtk_file(Grid*& gl, const int id, const std::vector<std::array<int,3>>& point_indices, std::string base_filename) {

    std::string filename = base_filename + "_" + std::to_string(id) + ".vtk";
    std::vector<std::array<double,3>> coords;
	std::array<double,3> icoords;
    coords.resize(point_indices.size());
    for(int pt=0;pt<point_indices.size();pt++){
        int i = point_indices[pt][0];
        int j = point_indices[pt][1];
        int k = point_indices[pt][2];
        icoords[0] = gl[id].x[i][j][k];
        icoords[1] = gl[id].y[i][j][k];
        icoords[2] = 1e-12;//gl[id].z[i][j][k];
		coords[pt] = icoords;
    }

    write_points_vtk(coords, filename.c_str());
}

void write_points_vtk(std::vector<std::array<double,3>>& coords, std::string filename)
{
        FILE* file_points_vtk;
        file_points_vtk = fopen(filename.c_str(),"w");
        fprintf(file_points_vtk, "%s\n","# vtk DataFile Version 3.0");
        fprintf(file_points_vtk, "%s\n","Wind turbine locations");
        fprintf(file_points_vtk, "%s\n","ASCII");
        fprintf(file_points_vtk, "%s\n","DATASET POLYDATA");
        fprintf(file_points_vtk, "%s %ld %s\n", "POINTS", coords.size(), "float");
        for(long unsigned int i=0; i<coords.size(); i++){
            fprintf(file_points_vtk, "%0.15g %0.15g %0.15g\n", coords[i][0], coords[i][1], coords[i][2]);
        }
        fclose(file_points_vtk);
}


