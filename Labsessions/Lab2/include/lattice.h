#pragma once
/* ****************************************************************************
 * Function that generates a fcc lattice in units of [Å]. Nc is the number of
 * primitive cells in each direction and a0 is the lattice parameter. The
 * positions of all the atoms are stored in pos which should be a matrix of the
 * size N x 3, where N is the number of atoms. The first, second and third
 * column correspond to the x,y and z coordinate respectively.
 *
 * Example usage:
 *
 * >>> init_fcc(pos, Nc, a0);
 *
 * ****************************************************************************/
void init_fcc(double **, int, double);


/* ****************************************************************************
 * Function that generates a sc lattice in units of [Å]. Nc is the number of
 * primitive cells in each direction and a0 is the lattice parameter. The
 * positions of all the atoms are stored in pos which should be a matrix of the
 * size N x 3, where N is the number of atoms. The first, second and third
 * column correspond to the x,y and z coordinate respectively.
 *
 * Example usage:
 *
 * >>> init_fcc(pos, Nc, a0);
 *
 * ****************************************************************************/
void init_sc(double **positions, int *material, int N, double lattice_param, unsigned int species);

void idx_nearest(double** idx_nearest, int sublattice, int n_atoms, int n_cells);
