#include <Grid.H>
#include <DataStruct.H>
using namespace std;

void Read_Grid(Grid*& gl, int& num_blocks, const std::string filename)
{
	double t;
	int Nblocks;

	//sprintf(filename, "full_hill_small_grid1.grd");
	ifstream myFile;
	myFile.open (filename.c_str(), ios::in | ios::binary);

	myFile.seekg(4, ios::cur);
	myFile.read((char*) &Nblocks, sizeof(int));
	myFile.seekg(4, ios::cur);
	cout << Nblocks << endl;

	gl = new Grid [Nblocks];
	num_blocks = Nblocks;
		
	myFile.seekg(4, ios::cur);
	for (int i=0;i<Nblocks;i++)
	{
		myFile.read((char*) &gl[i].Nx, sizeof(int));
		myFile.read((char*) &gl[i].Ny, sizeof(int));
		myFile.read((char*) &gl[i].Nz, sizeof(int));

		cout << gl[i].Nx << endl;
		cout << gl[i].Ny << endl;
		cout << gl[i].Nz << endl;

		gl[i].x = Create3DMatrix(gl[i].Nx,gl[i].Ny,gl[i].Nz);
		gl[i].y = Create3DMatrix(gl[i].Nx,gl[i].Ny,gl[i].Nz);
		gl[i].z = Create3DMatrix(gl[i].Nx,gl[i].Ny,gl[i].Nz);
		gl[i].iblank = Create3DMatrix_INT(gl[i].Nx,gl[i].Ny,gl[i].Nz);
	}
	myFile.seekg(4, ios::cur);

	for (int ii=0;ii<Nblocks;ii++)
	{
		myFile.seekg(4, ios::cur);
		for (int k=0;k<gl[ii].Nz;++k)
			for (int j=0;j<gl[ii].Ny;++j)
				for (int i=0;i<gl[ii].Nx;++i)
				{
					myFile.read((char*) &gl[ii].x[i][j][k], sizeof(t));
					//cout << gl[ii].x[i][j] << endl;
				}

		for (int k=0;k<gl[ii].Nz;++k)
			for (int j=0;j<gl[ii].Ny;++j)
				for (int i=0;i<gl[ii].Nx;++i)
				{
					myFile.read((char*) &gl[ii].y[i][j][k], sizeof(t));
					//cout << gl[ii].y[i][j] << endl;
				}

		for (int k=0;k<gl[ii].Nz;++k)
			for (int j=0;j<gl[ii].Ny;++j)
				for (int i=0;i<gl[ii].Nx;++i)
				{
					myFile.read((char*) &gl[ii].z[i][j][k], sizeof(t));
					//cout << "i = " << i << "    z = " << gl[ii].z[i][j][k] << endl;
				}

		for (int k=0;k<gl[ii].Nz;++k)
			for (int j=0;j<gl[ii].Ny;++j)
				for (int i=0;i<gl[ii].Nx;++i)
				{
					myFile.read((char*) &gl[ii].iblank[i][j][k], sizeof(int));
					//cout << gl[0].iblank[i][j] << endl;
				}

		myFile.seekg(4, ios::cur);

	}

	myFile.close();
	
	// WRITE THE GRID TO CHECK IF IT READ CORRECTLY
	/*ofstream outFile;
	outFile.open ("grid.xyz", ios::out | ios::binary);

	int size = 4;
	outFile.write((char*) &size, sizeof(int));
	outFile.write((char*) &Nblocks, sizeof(Nblocks));
	outFile.write((char*) &size, sizeof(int));

	size = 4*Nblocks*3;
	outFile.write((char*) &size, sizeof(int));
	for (int i=0;i<Nblocks;i++)
	{			
		outFile.write((char*) &gl[i].Nx, sizeof(gl[i].Nx));
		outFile.write((char*) &gl[i].Ny, sizeof(gl[i].Ny));
		outFile.write((char*) &gl[i].Nz, sizeof(gl[i].Nz));
	}
	outFile.write((char*) &size, sizeof(int));

	for (int ii=0;ii<Nblocks;ii++)
	{
		size = 8*3*gl[ii].Nx*gl[ii].Ny*gl[ii].Nz + 4*gl[ii].Nx*gl[ii].Ny*gl[ii].Nz;
		outFile.write((char*) &size, sizeof(int));
		for (int k=0;k<gl[ii].Nz;++k)
			for (int j=0;j<gl[ii].Ny;++j)
				for (int i=0;i<gl[ii].Nx;++i)
				{
					outFile.write((char*) &gl[ii].x[i][j][k], sizeof(gl[ii].x[i][j][k]));
				}

		for (int k=0;k<gl[ii].Nz;++k)
			for (int j=0;j<gl[ii].Ny;++j)
				for (int i=0;i<gl[ii].Nx;++i)
				{
					outFile.write((char*) &gl[ii].y[i][j][k], sizeof(gl[ii].y[i][j][k]));
				}

		for (int k=0;k<gl[ii].Nz;++k)
			for (int j=0;j<gl[ii].Ny;++j)
				for (int i=0;i<gl[ii].Nx;++i)
				{
					outFile.write((char*) &gl[ii].z[i][j][k], sizeof(gl[ii].z[i][j][k]));
				}

		for (int k=0;k<gl[ii].Nz;++k)
			for (int j=0;j<gl[ii].Ny;++j)
				for (int i=0;i<gl[ii].Nx;++i)
				{
					outFile.write((char*) &gl[ii].iblank[i][j][k], sizeof(int));
				}
		outFile.write((char*) &size, sizeof(int));
	}

	outFile.close();*/
}


void Deallocate_grid_memory(Grid* gl)
{
	int Nblocks = 2;
	for (int i=0;i<Nblocks;++i)
	{
		Delete3DMatrix(gl[i].x, gl[i].Nx, gl[i].Ny, gl[i].Nz);
		Delete3DMatrix(gl[i].y, gl[i].Nx, gl[i].Ny, gl[i].Nz);
		Delete3DMatrix(gl[i].z, gl[i].Nx, gl[i].Ny, gl[i].Nz);
		Delete3DMatrix_INT(gl[i].iblank, gl[i].Nx, gl[i].Ny, gl[i].Nz);
	}
}

