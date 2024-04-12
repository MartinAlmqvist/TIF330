#include <math.h>
void init_fcc(double **positions, int N, double lattice_param)
{
    double base[4][3] = {
	{0., 0., 0.},
	{0.5, 0.5, 0.},
	{0., 0.5, 0.5},
	{0.5, 0., 0.5}
    };
    int ind = 0;
    for (int x = 0; x < N; x++){
        for (int y = 0; y < N; y++){
            for (int z = 0; z < N; z++){
		for(int b = 0; b < 4; b++){
		    positions[ind][0] = lattice_param * (base[b][0] + x);
                    positions[ind][1] = lattice_param * (base[b][1] + y);
                    positions[ind][2] = lattice_param * (base[b][2] + z);
		    ind += 1;
		}
            }
        }
    }
}

void init_sc(double **positions, int *material, int N, double lattice_param, unsigned int species)
{
    double base[8][3] = {
	{0., 0., 0.},
	{1., 0., 0.},
	{0., 1., 0.},
	{1., 1., 0.},
	{0., 0., 1.},
	{1., 0., 1.},
	{0., 1., 1.},
	{1., 1., 1.},
    };
    int ind = 0;
    for (int x = 0; x < N; x++){
        for (int y = 0; y < N; y++){
            for (int z = 0; z < N; z++){
		for(int b = 0; b < 5; b++){
		    positions[ind][0] = lattice_param * (base[b][0] + x);
                    positions[ind][1] = lattice_param * (base[b][1] + y);
                    positions[ind][2] = lattice_param * (base[b][2] + z);
                    
                    // Set the type to 0 or 1 depending on the material
                    material[ind] = species;
                    
		    ind += 1;
		}
            }
        }
    }
}

void idx_nearest(double** idx_nearest, int sublattice, int n_atoms, int n_cells)
{
	// Nearest neighbour of atoms in A sublattice are all in B sublattice
	
	int n = pow(n_atoms, 1./3.);	// number of atoms in a row 
	int m = n * n; 		// number of atoms in a layer
	
	for (int a = 0; a <n_atoms; ++a){
		idx_nearest[a][0] = a;
		idx_nearest[a][1] = a - 1.;
		idx_nearest[a][2] = a - n;
		idx_nearest[a][3] = a - n - 1.;
		idx_nearest[a][4] = a - m;
		idx_nearest[a][5] = a - m - 1.;
		idx_nearest[a][6] = a - m - n;
		idx_nearest[a][7] = a - m - n - 1.;
		
		
		// Periodic boundary conditions for a sublattice A
		// x-direction
		if(a % n == 0){
			idx_nearest[a][1] += n;
			idx_nearest[a][2] += n;
			idx_nearest[a][5] += n;
			idx_nearest[a][7] += n;
		}
		// y direction
		if(a % m < n){
			idx_nearest[a][2] += n + 1;
			idx_nearest[a][3] += n + 1;
			idx_nearest[a][4] += n + 1;
			idx_nearest[a][5] += n + 1;
		}
		// z-direction: bottom of lattice to top of lattice
		if(a < m){
			idx_nearest[a][4] += m * n;
			idx_nearest[a][5] += m * n;
			idx_nearest[a][6] += m * n;
			idx_nearest[a][7] += m * n;
		}
		// z-direction B lattice
		//if(a > m-n == 0){
		//	idx_nearest[a][0] -= m * n;
		//}
		
		/*
		/for(int i = 0; i < 8; ++i){
			if (idx_nearest[a][i] < n_atoms){
				idx_nearest[a][i] += n_atoms;
			}
			if (idx_nearest[a][i] > n_atoms){
				idx_nearest[a][i] -= n_atoms;
			}
		}*/
		
	}	
}



