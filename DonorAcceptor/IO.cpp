#include <Grid.H>
#include <vector>
#include <iostream>
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


void write_points_vtk(std::vector<Point>& points, std::string filename)
{
        FILE* file_points_vtk;
        file_points_vtk = fopen(filename.c_str(),"w");
        fprintf(file_points_vtk, "%s\n","# vtk DataFile Version 3.0");
        fprintf(file_points_vtk, "%s\n","Wind turbine locations");
        fprintf(file_points_vtk, "%s\n","ASCII");
        fprintf(file_points_vtk, "%s\n","DATASET POLYDATA");
        fprintf(file_points_vtk, "%s %ld %s\n", "POINTS", points.size(), "float");
        for(long unsigned int i=0; i<points.size(); i++){
            fprintf(file_points_vtk, "%0.15g %0.15g %0.15g\n", points[i].x, points[i].y, points[i].z);
        }
        fclose(file_points_vtk);
}


