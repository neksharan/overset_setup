#include <Grid.H>
#include <vector>
#include <iostream>
using namespace std;
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


