#include <Grid.H>
#include <DataStruct.H>
using namespace std;

void Write_Grid(const double* xs, const double* ys, const double* zs,
				const double* xe, const double* ye, const double* ze,
				const int* Nx_pts, const int* Ny_pts, const int* Nz_pts,
				const int Nblocks, int base_grid_index, HoleCut& hole_cut, Grid* grid)
{
	for (int i=0;i<Nblocks;i++)
	{
		grid[i].Nx = Nx_pts[i];
		grid[i].Ny = Ny_pts[i];
		grid[i].Nz = Nz_pts[i];

		cout << grid[i].Nx << endl;
		cout << grid[i].Ny << endl;
		cout << grid[i].Nz << endl;

		grid[i].x = Create3DMatrix(grid[i].Nx,grid[i].Ny,grid[i].Nz);
		grid[i].y = Create3DMatrix(grid[i].Nx,grid[i].Ny,grid[i].Nz);
		grid[i].z = Create3DMatrix(grid[i].Nx,grid[i].Ny,grid[i].Nz);
		grid[i].iblank = Create3DMatrix_INT(grid[i].Nx,grid[i].Ny,grid[i].Nz);
	}

	for (int ii=0;ii<Nblocks;ii++)
	{
		if (ii==0)
		{
			double r[grid[ii].Nx];
			for (int i=0;i<grid[ii].Nx;++i)
			{	r[i] = xs[ii] + i * (xe[ii]-xs[ii])/(grid[ii].Nx-1);	}

			double theta[grid[ii].Ny];
			for (int j=0;j<grid[ii].Ny;++j)
			{	theta[j] = ys[ii] + j * (ye[ii]-ys[ii])/(grid[ii].Ny-1);	}

			for (int k=0;k<grid[ii].Nz;++k)
				for (int j=0;j<grid[ii].Ny;++j)
					for (int i=0;i<grid[ii].Nx;++i)
					{				
						grid[ii].x[i][j][k] = r[i]*cos(theta[j]);			
						grid[ii].y[i][j][k] = r[i]*sin(theta[j]);
						grid[ii].z[i][j][k] = 0.0;
						grid[ii].iblank[i][j][k] = 1;
					}
		}
		else
		{
			double dx = (xe[ii]-xs[ii])/(Nx_pts[ii]-1);
			double dy = (ye[ii]-ys[ii])/(Ny_pts[ii]-1);
			double dz = 0.0;
			for (int k=0;k<grid[ii].Nz;++k)
				for (int j=0;j<grid[ii].Ny;++j)
					for (int i=0;i<grid[ii].Nx;++i)
					{
						grid[ii].x[i][j][k] = xs[ii] + dx*i;	
						grid[ii].y[i][j][k] = ys[ii] + dy*j;		
						grid[ii].z[i][j][k] = zs[ii] + dz*k;
						grid[ii].iblank[i][j][k] = 1;
					}
		}
	}

	ofstream outFile;
	char filename[] = "cyl_acoustic_grid.xyz";
	outFile.open (filename, ios::out | ios::binary);

	int size = 4;
	outFile.write((char*) &size, sizeof(int));
	outFile.write((char*) &Nblocks, sizeof(Nblocks));
	outFile.write((char*) &size, sizeof(int));

	size = 4*Nblocks*3;
	outFile.write((char*) &size, sizeof(int));
	for (int i=0;i<Nblocks;i++)
	{			
		outFile.write((char*) &grid[i].Nx, sizeof(grid[i].Nx));
		outFile.write((char*) &grid[i].Ny, sizeof(grid[i].Ny));
		outFile.write((char*) &grid[i].Nz, sizeof(grid[i].Nz));
	}
	outFile.write((char*) &size, sizeof(int));

	for (int ii=0;ii<Nblocks;++ii)
	{
		size = 8*3*grid[ii].Nx*grid[ii].Ny*grid[ii].Nz + 4*grid[ii].Nx*grid[ii].Ny*grid[ii].Nz;
		outFile.write((char*) &size, sizeof(int));
		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
				{
					outFile.write((char*) &grid[ii].x[i][j][k], sizeof(grid[ii].x[i][j][k]));
				}

		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
				{
					outFile.write((char*) &grid[ii].y[i][j][k], sizeof(grid[ii].y[i][j][k]));
				}

		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
				{
					outFile.write((char*) &grid[ii].z[i][j][k], sizeof(grid[ii].z[i][j][k]));
				}

		for (int k=0;k<grid[ii].Nz;++k)
			for (int j=0;j<grid[ii].Ny;++j)
				for (int i=0;i<grid[ii].Nx;++i)
				{
					double xval = grid[ii].x[i][j][k];
					double yval = grid[ii].y[i][j][k];
					if (ii==base_grid_index && xval>=hole_cut.xmin && xval<=hole_cut.xmax 
						     && yval>=hole_cut.ymin && yval<=hole_cut.ymax)
					{	grid[ii].iblank[i][j][k] = 0;		}
					outFile.write((char*) &grid[ii].iblank[i][j][k], sizeof(int));
				}

		outFile.write((char*) &size, sizeof(int));
	}

	outFile.close();
	cout << "===============================" << endl;
	printf("Wrote \"%s\" \n", filename);
	cout << "===============================" << endl;
	Deallocate_grid_memory(grid);
}
