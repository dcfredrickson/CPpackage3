/*     nonlocal_byatom, part of the Fredrickson Group Chemical Pressure Package  

                  Copyright (C) 2012, by Daniel C. Fredrickson

                    Last modified:  Mar. 28, 2012

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>

#define CUT3D_COMMAND "cut3d"
#define PI 3.14159265
#define NEQVOX 20
#define NIONMAX 200
#define TWO_PI 6.28318531
#define ONE_TWO_PI 0.159154943
#define PI_5_4 4.182513398
#define MAX_YAEHMOPFILES 10
#define LONGEST_FILENAME 100
#define NATOMS_MAX 300

#define TBANDS_MAX 3000
#define BANDS_MAX 1000
#define NORBS_MAX 3000
#define MAX_SYM 1000
#define KPOINTS_MAX 1000
#define WAVES_MAX 301000
#define NTYPES_MAX 10
#define R_MAX 10.0 /* bohr */
#define R_MULT 6.0
#define R_BOHR 0.52917721092 /* angstrom */

int ngx, ngy, ngz;
int kam, kbm, kcm;
int nsym;
int symrel[3][3][MAX_SYM];
double tnons[3][MAX_SYM];
int l_max = 0;

struct psp_file {
    double rloc[NTYPES_MAX];
    double rrs[NTYPES_MAX];
    double rrp[NTYPES_MAX];
    double rrd[NTYPES_MAX];
    double rrf[NTYPES_MAX];
    double cc1[NTYPES_MAX];
    double cc2[NTYPES_MAX];
    double cc3[NTYPES_MAX];
    double cc4[NTYPES_MAX];
    double h[3][3][4][NTYPES_MAX];
    double h11s[NTYPES_MAX];
    double h22s[NTYPES_MAX];
    double h33s[NTYPES_MAX];
    double h11p[NTYPES_MAX];
    double h22p[NTYPES_MAX];
    double h33p[NTYPES_MAX];
    double h11d[NTYPES_MAX];
    double h22d[NTYPES_MAX];
    double h33d[NTYPES_MAX];
    double h11f[NTYPES_MAX];
    double h22f[NTYPES_MAX];
    double h33f[NTYPES_MAX];
    double k11p[NTYPES_MAX];
    double k22p[NTYPES_MAX];
    double k33p[NTYPES_MAX];
    double k11d[NTYPES_MAX];
    double k22d[NTYPES_MAX];
    double k33d[NTYPES_MAX];
    double k11f[NTYPES_MAX];
    double k22f[NTYPES_MAX];
    double k33f[NTYPES_MAX];
    double k12p[NTYPES_MAX];
    double k23p[NTYPES_MAX];
    double k13p[NTYPES_MAX];
    double k12d[NTYPES_MAX];
    double k23d[NTYPES_MAX];
    double k13d[NTYPES_MAX];
    double k12f[NTYPES_MAX];
    double k23f[NTYPES_MAX];
    double k13f[NTYPES_MAX];
    int l_max[NTYPES_MAX];
    double r_max[NTYPES_MAX];
}
psp;

struct wfk_file {
    int nbands;
    int nkpts;
    int npw[2][KPOINTS_MAX];
    int nspinor[2][KPOINTS_MAX];
    int nband[2][KPOINTS_MAX];
    int kg[WAVES_MAX][3];
    double eigen[BANDS_MAX];
    double occ[BANDS_MAX];
    double cg[WAVES_MAX][BANDS_MAX][2];
}
wfk1;

double nonlocalE_byatom[3][NATOMS_MAX][6];

struct XSFfile {
    char systemname[100];
    double cellvolume;
    double scale_xyz;
    double cella_x, cella_y, cella_z;
    double cellb_x, cellb_y, cellb_z;
    double cellc_x, cellc_y, cellc_z;
    int atomicno[NATOMS_MAX];
    double Xcart[NATOMS_MAX];
    double Ycart[NATOMS_MAX];
    double Zcart[NATOMS_MAX];
    int NIONS;
    int NGX, NGY, NGZ;
    double ** * grid;
    double VoxelV;
    double NELECT;
    double epsatm[NATOMS_MAX];
    int epsatm_map[NATOMS_MAX];
    int epsatm_map_type[NATOMS_MAX];
    double nelectrons;
}
den0, den1, nlden1, nlden2;

double cgt[2][WAVES_MAX][2];
double ga[WAVES_MAX], gb[WAVES_MAX], gc[WAVES_MAX];
double gx[WAVES_MAX], gy[WAVES_MAX], gz[WAVES_MAX];
double costheta[WAVES_MAX];
double theta[WAVES_MAX];
double phi[WAVES_MAX];
double g[WAVES_MAX];
double Plm[WAVES_MAX][4][4];
/* p1=projp(l,i,g[pw1],pspin,atom_type,cellvolume); */
double p[WAVES_MAX][4][3][NTYPES_MAX];
int typat[NATOMS_MAX];

struct ContactVol {
    int ** * neighcount2, ** * ionmap[NEQVOX];
    double ** * swj;
    double ** * Plj_RE[NEQVOX][4][8][4], ** * Plj_IM[NEQVOX][4][8][3];
}
vmap;

FILE * inputfile;
FILE * outputfile;
char filename[LONGEST_FILENAME];
int k_use = -5;

int AllocDbl(struct XSFfile * gridin) {
    /* called by: */
    /* calls: none */
    int gridx;
    int gridy;
    int gridz;
    int jx, jy, jz;
    gridx = gridin -> NGX;
    gridy = gridin -> NGY;
    gridz = gridin -> NGZ;
    gridin -> grid = (double ** * ) malloc(gridx * sizeof(double ** ));
    for (jx = 0; jx < gridx; jx++) {
        gridin -> grid[jx] = (double ** ) malloc(gridy * sizeof(double * ));
        for (jy = 0; jy < gridy; jy++) {
            gridin -> grid[jx][jy] = (double * ) malloc(gridz * sizeof(double));
        }
    }
    return 0;
}

int AllocInt(struct ContactVol * map) {
    /* called by: main */
    /* calls: none */
    int gridx = ngx + 1, gridy = ngy + 1, gridz = ngz + 1, i = 0, jx = 0, jy = 0, m, l, j;
    printf("  Allocating %d x %d x %d grids for projectors.\n", ngx, ngy, ngz);
    map -> neighcount2 = (int ** * ) malloc(gridx * sizeof(int ** ));
    map -> swj = (double ** * ) malloc(gridx * sizeof(double ** ));
    for (jx = 0; jx < gridx; jx++) {
        map -> neighcount2[jx] = (int ** ) malloc(gridy * sizeof(int * ));
        map -> swj[jx] = (double ** ) malloc(gridy * sizeof(double * ));
        for (jy = 0; jy < gridy; jy++) {
            map -> neighcount2[jx][jy] = (int * ) malloc(gridz * sizeof(int));
            map -> swj[jx][jy] = (double * ) malloc(gridz * sizeof(double));
        }
    }
    for (i = 0; i < NEQVOX; i++) {
        map -> ionmap[i] = (int ** * ) malloc(gridx * sizeof(int ** ));
        for (jx = 0; jx < gridx; jx++) {
            map -> ionmap[i][jx] = (int ** ) malloc(gridy * sizeof(int * ));
            for (jy = 0; jy < gridy; jy++) {
                map -> ionmap[i][jx][jy] = (int * ) malloc(gridz * sizeof(int));
            }
        }
    }
    for (i = 0; i < NEQVOX; i++) {
        for (l = 0; l < l_max; l++) {
            for (m = 0; m < 2 * (l_max - 1) + 2; m++) {
                for (j = 0; j < 3; j++) {
                    map -> Plj_RE[i][l][m][j] = (double ** * ) malloc(gridx * sizeof(double ** ));
                    map -> Plj_IM[i][l][m][j] = (double ** * ) malloc(gridx * sizeof(double ** ));
                    for (jx = 0; jx < gridx; jx++) {
                        map -> Plj_RE[i][l][m][j][jx] = (double ** ) malloc(gridy * sizeof(double * ));
                        map -> Plj_IM[i][l][m][j][jx] = (double ** ) malloc(gridy * sizeof(double * ));
                        for (jy = 0; jy < gridy; jy++) {
                            map -> Plj_RE[i][l][m][j][jx][jy] = (double * ) malloc(gridz * sizeof(double));
                            map -> Plj_IM[i][l][m][j][jx][jy] = (double * ) malloc(gridz * sizeof(double));
                        }
                    }
                }
            }
        }
    }
    return 0;
}

int Getkm(struct XSFfile * gridin) {
    /* called by main */
    /* calls: none */
    double fa = 0.0, fb = 0.0, fc = 0.0;
    fa = pow((gridin -> cellb_y * gridin -> cellc_z - gridin -> cellb_z * gridin -> cellc_y), 2) +
        pow((gridin -> cellb_z * gridin -> cellc_x - gridin -> cellb_x * gridin -> cellc_z), 2) +
        pow((gridin -> cellb_x * gridin -> cellc_y - gridin -> cellb_y * gridin -> cellc_x), 2);
    fb = pow((gridin -> cellc_y * gridin -> cella_z - gridin -> cellc_z * gridin -> cella_y), 2) +
        pow((gridin -> cellc_z * gridin -> cella_x - gridin -> cellc_x * gridin -> cella_z), 2) +
        pow((gridin -> cellc_x * gridin -> cella_y - gridin -> cellc_y * gridin -> cella_x), 2);
    fc = pow((gridin -> cella_y * gridin -> cellb_z - gridin -> cella_z * gridin -> cellb_y), 2) +
        pow((gridin -> cella_z * gridin -> cellb_x - gridin -> cella_x * gridin -> cellb_z), 2) +
        pow((gridin -> cella_x * gridin -> cellb_y - gridin -> cella_y * gridin -> cellb_x), 2);
    kam = (int) ceil(R_MAX / (gridin -> cellvolume / sqrt(fa)));
    kbm = (int) ceil(R_MAX / (gridin -> cellvolume / sqrt(fb)));
    kcm = (int) ceil(R_MAX / (gridin -> cellvolume / sqrt(fc)));
    printf("  Using supercell range %d %d %d\n", kam, kbm, kcm);
    if (kam > 3 || kbm > 3 || kcm > 3) {
        printf("\n  BAD NEWS: Unit cell too small or distorted!\n");
        return 1;
    }
    return 0;
}

void xyz2sph(double X, double Y, double Z, double * r, double * theta, double * phi) {
    /* Z = r cos(theta), X = r sin(theta) cos(phi) Y = r sin(theta)sin(phi) */
    * r = pow((pow(X, 2.0) + pow(Y, 2.0) + pow(Z, 2.0)), 0.5);
    if ( * r > 0.0) {
        * theta = acos(Z / * r);
        if (fabs(sin( * theta)) > 0.0) {
            * phi = acos(X / ( * r * sin( * theta)));
            if (X / ( * r * sin( * theta)) > 1.0) {
                * phi = acos(1.0);
            }
            if (X / ( * r * sin( * theta)) < -1.0) {
                * phi = acos(-1.0);
            }
            if (Y < 0.0) * phi = - * phi;
        } else * phi = 0;
    } else {
        * phi = 0.0;
        * theta = 0.0;
    }
}

double projp_r(int l, int i, double r, struct psp_file * pspin, int typat) {
    double p = 0.0;
    double r_l;
    if (l == 0) r_l = pspin -> rrs[typat];
    if (l == 1) r_l = pspin -> rrp[typat];
    if (l == 2) r_l = pspin -> rrd[typat];
    if (l == 3) r_l = pspin -> rrf[typat];
    if (r_l > 0.0) {
        p = sqrt(2.0) * pow(r, 1.0 * l + 2.0 * (1.0 * i + 1.0) - 2.0) * exp(-pow(r, 2.0) / (2.0 * pow(r_l, 2.0))) / (pow(r_l, 1.0 * l + 2.0 * (1.0 * i + 1.0) - 0.5) * sqrt(gsl_sf_gamma(l + 2.0 * (1.0 * i + 1.0) - 0.5)));
    }
    return (p);
}

int CoordSearchSphere(struct XSFfile * gridin, struct ContactVol * map) {
    /* called by: main */
    /* calls: Getwj */
    int atom = 0, atom2 = 0, atom3[NEQVOX], index = 0, index2, i = 0, j = 0, k, jx = 0, jy = 0, jz = 0, l, m, atom_type;
    int ka1 = 0, kb1 = 0, kc1 = 0, ka2 = 0, kb2 = 0, kc2 = 0, ka3[NEQVOX], kb3[NEQVOX], kc3[NEQVOX];
    int ka = 0, kb = 0, kc = 0, check[NEQVOX], ngcount = 0, ngp = 0, ngp0 = 0, voxtot = ngx * ngy * ngz;
    double interatom_dist, dist = 0.0, distmin, dmax = 0.0, voxcenter_x = 0.0, voxcenter_y = 0.0, voxcenter_z = 0.0;
    double wj_temp[NEQVOX], wmax = 0.0, wmax2 = 0.0, xf = 0.0, yf = 0.0, zf = 0.0, xc = 0.0, yc = 0.0, zc = 0.0, xc2 = 0.0, yc2 = 0.0, zc2 = 0.0;
    double cella = 0.0, cellb = 0.0, cellc = 0.0, minstep = 0.0, tolerance = 0.0, r_max;
    double stepx, stepy, stepz, testx, testy, testz, gradx, grady, gradz, grad_norm;
    double atom_volumes[NIONMAX], atom_charges[NIONMAX];
    double xf2, yf2, zf2;
    int stop, step, jx2, jy2, jz2, jx2c, jy2c, jz2c;
    double weightmatrix[NEQVOX][NEQVOX];
    double weightmatrix2[NEQVOX][NEQVOX];
    int neighborindex[NEQVOX];
    double weightmatrix_max;
    int contactcounter, counter;
    int stepmax;
    int step_upper, step_lower, step_new;
    int ionmap_temp, ionmap_temp2;
    int ionmap_list[NIONMAX * 27];
    double ionmap_dist[NIONMAX * 27], diffx, diffy, diffz;
    char bader_name[100];
    int foundit;
    int shared_voxels, symm_op;
    int ka0, kb0, kc0, atom0;
    double sinmphi[4], cosmphi[4], costheta1, phi, theta;
    cella = pow(gridin -> cella_x * gridin -> cella_x + gridin -> cella_y * gridin -> cella_y + gridin -> cella_z * gridin -> cella_z, 0.5) / (ngx);
    cellb = pow(gridin -> cellb_x * gridin -> cellb_x + gridin -> cellb_y * gridin -> cellb_y + gridin -> cellb_z * gridin -> cellb_z, 0.5) / (ngy);
    cellc = pow(gridin -> cellc_x * gridin -> cellc_x + gridin -> cellc_y * gridin -> cellc_y + gridin -> cellc_z * gridin -> cellc_z, 0.5) / (ngz);
    minstep = cella;
    if (cellb < minstep) minstep = cellb;
    if (cellc < minstep) minstep = cellc;
    tolerance = 0.0004 * minstep;
    /* every voxel in the unit cell */
    for (jz = 0; jz < ngz; jz++) {
        for (jy = 0; jy < ngy; jy++) {
            for (jx = 0; jx < ngx; jx++) {
                map -> neighcount2[jx][jy][jz] = 0;
                map -> swj[jx][jy][jz] = 0.0;
            }
        }
    }
    printf("Creating real space projectors...\n");
    for (atom = 0; atom < gridin -> NIONS; atom++) {
        atom_volumes[atom] = 0.0;
        atom_charges[atom] = 0.0;
    }
    for (jz = 0; jz < ngz; jz++) {
        /* fractional coordinates */
        zf = (double) jz / (double) ngz;
        for (jy = 0; jy < ngy; jy++) {
            yf = (double) jy / (double) ngy;
            for (jx = 0; jx < ngx; jx++) {
                xf = (double) jx / (double) ngx;
                voxcenter_x = xf * gridin -> cella_x + yf * gridin -> cellb_x + zf * gridin -> cellc_x;
                voxcenter_y = xf * gridin -> cella_y + yf * gridin -> cellb_y + zf * gridin -> cellc_y;
                voxcenter_z = xf * gridin -> cella_z + yf * gridin -> cellb_z + zf * gridin -> cellc_z;
                for (atom = 0; atom < gridin -> NIONS; atom++) {
                    atom_type = typat[atom] - 1;
                    for (ka = -kam; ka <= kam; ka++) {
                        for (kb = -kbm; kb <= kbm; kb++) {
                            for (kc = -kcm; kc <= kcm; kc++) {
                                xc = gridin -> Xcart[atom] + ka * gridin -> cella_x + kb * gridin -> cellb_x + kc * gridin -> cellc_x;
                                yc = gridin -> Ycart[atom] + ka * gridin -> cella_y + kb * gridin -> cellb_y + kc * gridin -> cellc_y;
                                zc = gridin -> Zcart[atom] + ka * gridin -> cella_z + kb * gridin -> cellb_z + kc * gridin -> cellc_z;
                                dist = sqrt((voxcenter_x - xc) * (voxcenter_x - xc) + (voxcenter_y - yc) * (voxcenter_y - yc) + (voxcenter_z - zc) * (voxcenter_z - zc));
                                if (dist <= R_MULT * psp.r_max[atom_type]) {
                                    index = map -> neighcount2[jx][jy][jz];
                                    if (index == NEQVOX) {
                                        printf("Voxel has too many neighbors.  Increase NEQVOX or decrease R_MULT.\n");
                                        exit(0);
                                    }
                                    map -> ionmap[index][jx][jy][jz] = ((ka + 3) << 13) + ((kb + 3) << 10) + ((kc + 3) << 7) + atom;
                                    map -> swj[jx][jy][jz] += 1.0;
                                    map -> neighcount2[jx][jy][jz]++;
                                    if (index > 0) shared_voxels++;
                                    atom_type = typat[atom] - 1;
                                    xyz2sph(voxcenter_x - xc, voxcenter_y - yc, voxcenter_z - zc, & dist, & theta, & phi);
                                    costheta1 = cos(theta);
                                    for (m = 1; m < 4; m++) {
                                        sinmphi[m] = sin(m * phi);
                                        cosmphi[m] = cos(m * phi);
                                    }
                                    for (j = 0; j < 3; j++) {
                                        map -> Plj_RE[index][0][0][j][jx][jy][jz] = projp_r(0, j, dist, & psp, atom_type) * gsl_sf_legendre_sphPlm(0, 0, 0);
                                        map -> Plj_IM[index][0][0][j][jx][jy][jz] = 0.0;
                                        for (l = 1; l < l_max; l++) {
                                            m = 0;
                                            map -> Plj_RE[index][l][m][j][jx][jy][jz] = projp_r(l, j, dist, & psp, atom_type) * gsl_sf_legendre_sphPlm(l, m, costheta1);
                                            map -> Plj_IM[index][l][m][j][jx][jy][jz] = projp_r(l, j, dist, & psp, atom_type) * gsl_sf_legendre_sphPlm(l, m, costheta1);
                                            for (m = 1; m < l + 1; m++) {
                                                map -> Plj_RE[index][l][2 * m + 1][j][jx][jy][jz] = projp_r(l, j, dist, & psp, atom_type) * gsl_sf_legendre_sphPlm(l, m, costheta1) * cosmphi[m];
                                                map -> Plj_RE[index][l][2 * m][j][jx][jy][jz] = projp_r(l, j, dist, & psp, atom_type) * gsl_sf_legendre_sphPlm(l, m, costheta1) * cosmphi[m];
                                                map -> Plj_IM[index][l][2 * m + 1][j][jx][jy][jz] = projp_r(l, j, dist, & psp, atom_type) * gsl_sf_legendre_sphPlm(l, m, costheta1) * sinmphi[m];
                                                map -> Plj_IM[index][l][2 * m][j][jx][jy][jz] = -projp_r(l, j, dist, & psp, atom_type) * gsl_sf_legendre_sphPlm(l, m, costheta1) * sinmphi[m];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                for (i = 0; i < map -> neighcount2[jx][jy][jz]; i++) {
                    atom = map -> ionmap[i][jx][jy][jz] & 127;
                    atom_volumes[atom] += gridin -> VoxelV / map -> swj[jx][jy][jz];
                    atom_charges[atom] += gridin -> grid[jx][jy][jz] * gridin -> VoxelV / map -> swj[jx][jy][jz];
                }
            }
        }
    }
    ngcount = 0;
    printf("\n  Sphere integration results \n");
    for (atom = 0; atom < gridin -> NIONS; atom++) {
        printf("   Atom %d:  charge:  %8.5lf   volume:  %8.5lf  \n", atom + 1, atom_charges[atom], atom_volumes[atom]);
    }
    printf("\n ...Finished\n");
    return 0;
}

int CoordSearchSphere2(struct XSFfile * gridin, struct ContactVol * map) {
    /* called by: main */
    /* calls: Getwj */
    int atom = 0, atom2 = 0, atom3[NEQVOX], index = 0, index2, i = 0, j = 0, k, jx = 0, jy = 0, jz = 0, l, m, atom_type;
    int ka1 = 0, kb1 = 0, kc1 = 0, ka2 = 0, kb2 = 0, kc2 = 0, ka3[NEQVOX], kb3[NEQVOX], kc3[NEQVOX];
    int ka = 0, kb = 0, kc = 0, check[NEQVOX], ngcount = 0, ngp = 0, ngp0 = 0, voxtot = ngx * ngy * ngz;
    double interatom_dist, dist = 0.0, distmin, dmax = 0.0, voxcenter_x = 0.0, voxcenter_y = 0.0, voxcenter_z = 0.0;
    double wj_temp[NEQVOX], wmax = 0.0, wmax2 = 0.0, xf = 0.0, yf = 0.0, zf = 0.0, xc = 0.0, yc = 0.0, zc = 0.0, xc2 = 0.0, yc2 = 0.0, zc2 = 0.0;
    double cella = 0.0, cellb = 0.0, cellc = 0.0, minstep = 0.0, tolerance = 0.0, r_max;
    double stepx, stepy, stepz, testx, testy, testz, gradx, grady, gradz, grad_norm;
    double atom_volumes[NIONMAX], atom_charges[NIONMAX];
    double xf2, yf2, zf2;
    int stop, step, jx2, jy2, jz2, jx2c, jy2c, jz2c;
    double weightmatrix[NEQVOX][NEQVOX];
    double weightmatrix2[NEQVOX][NEQVOX];
    int neighborindex[NEQVOX];
    double weightmatrix_max;
    int contactcounter, counter;
    int stepmax;
    int step_upper, step_lower, step_new;
    int ionmap_temp, ionmap_temp2;
    int ionmap_list[NIONMAX * 27];
    double ionmap_dist[NIONMAX * 27], diffx, diffy, diffz;
    char bader_name[100];
    int foundit;
    int shared_voxels, symm_op;
    int ka0, kb0, kc0, atom0;
    double sinmphi[4], cosmphi[4], costheta1, phi, theta;
    cella = pow(gridin -> cella_x * gridin -> cella_x + gridin -> cella_y * gridin -> cella_y + gridin -> cella_z * gridin -> cella_z, 0.5) / (ngx);
    cellb = pow(gridin -> cellb_x * gridin -> cellb_x + gridin -> cellb_y * gridin -> cellb_y + gridin -> cellb_z * gridin -> cellb_z, 0.5) / (ngy);
    cellc = pow(gridin -> cellc_x * gridin -> cellc_x + gridin -> cellc_y * gridin -> cellc_y + gridin -> cellc_z * gridin -> cellc_z, 0.5) / (ngz);
    minstep = cella;
    if (cellb < minstep) minstep = cellb;
    if (cellc < minstep) minstep = cellc;
    tolerance = 0.0004 * minstep;
    /* every voxel in the unit cell */
    for (jz = 0; jz < ngz; jz++) {
        for (jy = 0; jy < ngy; jy++) {
            for (jx = 0; jx < ngx; jx++) {
                map -> neighcount2[jx][jy][jz] = 0;
                map -> swj[jx][jy][jz] = 0.0;
            }
        }
    }
    printf("Creating real space projectors...\n");
    for (atom = 0; atom < gridin -> NIONS; atom++) {
        atom_volumes[atom] = 0.0;
        atom_charges[atom] = 0.0;
    }
    for (jz = 0; jz < ngz; jz++) {
        /* fractional coordinates */
        zf = (double) jz / (double) ngz;
        for (jy = 0; jy < ngy; jy++) {
            yf = (double) jy / (double) ngy;
            for (jx = 0; jx < ngx; jx++) {
                xf = (double) jx / (double) ngx;
                voxcenter_x = xf * gridin -> cella_x + yf * gridin -> cellb_x + zf * gridin -> cellc_x;
                voxcenter_y = xf * gridin -> cella_y + yf * gridin -> cellb_y + zf * gridin -> cellc_y;
                voxcenter_z = xf * gridin -> cella_z + yf * gridin -> cellb_z + zf * gridin -> cellc_z;
                for (atom = 0; atom < gridin -> NIONS; atom++) {
                    atom_type = typat[atom] - 1;
                    for (ka = -kam; ka <= kam; ka++) {
                        for (kb = -kbm; kb <= kbm; kb++) {
                            for (kc = -kcm; kc <= kcm; kc++) {
                                xc = gridin -> Xcart[atom] + ka * gridin -> cella_x + kb * gridin -> cellb_x + kc * gridin -> cellc_x;
                                yc = gridin -> Ycart[atom] + ka * gridin -> cella_y + kb * gridin -> cellb_y + kc * gridin -> cellc_y;
                                zc = gridin -> Zcart[atom] + ka * gridin -> cella_z + kb * gridin -> cellb_z + kc * gridin -> cellc_z;
                                dist = sqrt((voxcenter_x - xc) * (voxcenter_x - xc) + (voxcenter_y - yc) * (voxcenter_y - yc) + (voxcenter_z - zc) * (voxcenter_z - zc));
                                if (dist <= R_MULT * psp.r_max[atom_type]) {
                                    index = map -> neighcount2[jx][jy][jz];
                                    if (index == NEQVOX) {
                                        printf("Voxel has too many neighbors.  Increase NEQVOX or decrease R_MULT.\n");
                                        exit(0);
                                    }
                                    map -> ionmap[index][jx][jy][jz] = ((ka + 3) << 13) + ((kb + 3) << 10) + ((kc + 3) << 7) + atom;
                                    map -> swj[jx][jy][jz] += 1.0;
                                    map -> neighcount2[jx][jy][jz]++;
                                    if (index > 0) shared_voxels++;
                                    for (j = 0; j < 3; j++) {
                                        map -> Plj_RE[index][0][0][j][jx][jy][jz] = projp_r(0, j, dist, & psp, atom_type);
                                        map -> Plj_IM[index][0][0][j][jx][jy][jz] = 0.0;
                                        for (l = 1; l < l_max; l++) {
                                            m = 0;
                                            map -> Plj_RE[index][l][m][j][jx][jy][jz] = projp_r(l, j, dist, & psp, atom_type);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                for (i = 0; i < map -> neighcount2[jx][jy][jz]; i++) {
                    atom = map -> ionmap[i][jx][jy][jz] & 127;
                    atom_volumes[atom] += gridin -> VoxelV / map -> swj[jx][jy][jz];
                    atom_charges[atom] += gridin -> grid[jx][jy][jz] * gridin -> VoxelV / map -> swj[jx][jy][jz];
                }
            }
        }
    }
    ngcount = 0;
    printf("\n  Sphere integration results \n");
    for (atom = 0; atom < gridin -> NIONS; atom++) {
        printf("   Atom %d:  charge:  %8.5lf   volume:  %8.5lf  \n", atom + 1, atom_charges[atom], atom_volumes[atom]);
    }
    printf("\n ...Finished\n");
    return 0;
}

void finish_line(FILE * f3) {
    int cont = 0;
    char check;
    while (cont == 0) {
        check = getc(f3);
        if ((check == 10) || (check == EOF)) cont = 1;
    }
}

void outputXSF(struct XSFfile * XSFIN, struct XSFfile * XSFOUT, FILE * f2) {
    int jx, jy, jz;
    int j;
    int line_counter;
    fprintf(f2, " DIM-GROUP\n");
    fprintf(f2, " 3  1\n");
    fprintf(f2, " PRIMVEC\n");
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cella_x * R_BOHR, XSFIN -> cella_y * R_BOHR, XSFIN -> cella_z * R_BOHR);
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cellb_x * R_BOHR, XSFIN -> cellb_y * R_BOHR, XSFIN -> cellb_z * R_BOHR);
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cellc_x * R_BOHR, XSFIN -> cellc_y * R_BOHR, XSFIN -> cellc_z * R_BOHR);
    fprintf(f2, " PRIMCOORD\n");
    fprintf(f2, "%12d%3d\n", XSFIN -> NIONS, 1);
    for (j = 0; j < XSFIN -> NIONS; j++) {
        fprintf(f2, "%9d%20.10lf%20.10lf%20.10lf\n", XSFIN -> atomicno[j], XSFIN -> Xcart[j] * R_BOHR, XSFIN -> Ycart[j] * R_BOHR, XSFIN -> Zcart[j] * R_BOHR);
    }
    fprintf(f2, " ATOMS\n");
    for (j = 0; j < XSFIN -> NIONS; j++) {
        fprintf(f2, "%9d%20.10lf%20.10lf%20.10lf\n", XSFIN -> atomicno[j], XSFIN -> Xcart[j] * R_BOHR, XSFIN -> Ycart[j] * R_BOHR, XSFIN -> Zcart[j] * R_BOHR);
    }
    fprintf(f2, " BEGIN_BLOCK_DATAGRID3D\n");
    fprintf(f2, " datagrids\n");
    fprintf(f2, " DATAGRID_3D_DENSITY\n");
    fprintf(f2, "%12d%12d%12d\n", XSFIN -> NGX, XSFIN -> NGY, XSFIN -> NGZ);
    fprintf(f2, " 0.0 0.0 0.0\n");
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cella_x * R_BOHR, XSFIN -> cella_y * R_BOHR, XSFIN -> cella_z * R_BOHR);
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cellb_x * R_BOHR, XSFIN -> cellb_y * R_BOHR, XSFIN -> cellb_z * R_BOHR);
    fprintf(f2, "%17.10lf%19.10lf%19.10lf\n", XSFIN -> cellc_x * R_BOHR, XSFIN -> cellc_y * R_BOHR, XSFIN -> cellc_z * R_BOHR);
    line_counter = 0;
    for (jz = 0; jz < XSFIN -> NGZ; jz++) {
        for (jy = 0; jy < XSFIN -> NGY; jy++) {
            for (jx = 0; jx < XSFIN -> NGX; jx++) {
                line_counter++;
                fprintf(f2, "%20.10lf", XSFOUT -> grid[jx][jy][jz]);
                if (line_counter == 6) {
                    fprintf(f2, "\n");
                    line_counter = 0;
                }
            }
        }
    }
    fprintf(f2, " END_DATAGRID_3D\n");
    fprintf(f2, " END_BLOCK_DATAGRID3D\n");

}

void read_psp_data(char filename[200], struct psp_file * pspin) {
    FILE * f2;
    int jx, jy, jz, ha, hb, hc;
    int j, l;
    int line_counter;
    int nvoxels_atm;
    int check;
    int stop = 0;
    int stop2;
    char str[1000];
    int NGX, NGY, NGZ, k, ixc = -10;
    double cell_volume, dist;
    double rcut;
    double xf, yf, zf;
    double voxel_centerx;
    double voxel_centery;
    double voxel_centerz;
    int n_epsatm_values = 0;
    double Z;
    int typat;
    double Z_ion_temp;
    double Z_ion[NATOMS_MAX];
    int Z_values[NATOMS_MAX];
    double epsatm_values[NATOMS_MAX];
    ixc = 1;
    f2 = fopen(filename, "r");
    typat = 0;
    while (stop == 0) {
        check = fscanf(f2, "%s", str);
        if (check == EOF) {
            stop = 1;
        }
        if (strcmp(str, "ixc") == 0) {
            fscanf(f2, "%d", & ixc);
            printf("XC correlation functional = %d\n", ixc);
        }
        if (strcmp(str, "ETOT") == 0) stop = 1;
        if ((strcmp(str, "pspini:") == 0) && (ixc == 1)) {
            finish_line(f2);
            finish_line(f2);
            finish_line(f2);
            fscanf(f2, "%s", str);
            fscanf(f2, "%lf %lf", & Z, & Z_ion[typat]);
            Z_values[typat] = Z;
            pspin -> r_max[typat] = 0.0;
            stop2 = 0;
            while (stop2 == 0) {
                check = fscanf(f2, "%s", str);

                if (strcmp(str, "rloc=") == 0) {
                    fscanf(f2, "%lf", & pspin -> rloc[typat]);
                }
                if (strcmp(str, "cc1") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc1[typat]);
                }
                if (strcmp(str, "cc2") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc2[typat]);
                }
                if (strcmp(str, "cc3") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc3[typat]);
                }
                if (strcmp(str, "cc4") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc4[typat]);
                }
                if (strcmp(str, "rrs") == 0) {
                    fscanf(f2, " = %lf", & pspin -> rrs[typat]);
                    if (pspin -> rrs[typat] > 0) pspin -> l_max[typat] = 1;
                    if (pspin -> rrs[typat] > pspin -> r_max[typat]) pspin -> r_max[typat] = pspin -> rrs[typat];
                }
                if (strcmp(str, "h11s=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[0][0][0][typat]);
                }
                if (strcmp(str, "h22s=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[1][1][0][typat]);
                }
                if (strcmp(str, "h33s=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[2][2][0][typat]);
                }
                if (strcmp(str, "rrp") == 0) {
                    fscanf(f2, " = %lf", & pspin -> rrp[typat]);
                    if (pspin -> rrp[typat] > 0) pspin -> l_max[typat] = 2;
                    if (pspin -> rrp[typat] > pspin -> r_max[typat]) pspin -> r_max[typat] = pspin -> rrp[typat];
                }
                if (strcmp(str, "h11p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[0][0][1][typat]);
                }
                if (strcmp(str, "h22p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[1][1][1][typat]);
                }
                if (strcmp(str, "h33p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[2][2][1][typat]);
                }
                if (strcmp(str, "k11p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k11p[typat]);
                }
                if (strcmp(str, "k22p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k22p[typat]);
                }
                if (strcmp(str, "k33p=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k33p[typat]);
                }
                if (strcmp(str, "rrd") == 0) {
                    fscanf(f2, " = %lf", & pspin -> rrd[typat]);
                    if (pspin -> rrd[typat] > 0) pspin -> l_max[typat] = 3;
                    if (pspin -> rrd[typat] > pspin -> r_max[typat]) pspin -> r_max[typat] = pspin -> rrd[typat];
                }
                if (strcmp(str, "h11d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[0][0][2][typat]);
                }
                if (strcmp(str, "h22d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[1][1][2][typat]);
                }
                if (strcmp(str, "h33d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[2][2][2][typat]);
                }
                if (strcmp(str, "k11d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k11d[typat]);
                }
                if (strcmp(str, "k22d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k22d[typat]);
                }
                if (strcmp(str, "k33d=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k33d[typat]);
                }
                if (strcmp(str, "rrf") == 0) {
                    fscanf(f2, " = %lf", & pspin -> rrf[typat]);
                    if (pspin -> rrf[typat] > 0) pspin -> l_max[typat] = 4;
                    if (pspin -> rrf[typat] > pspin -> r_max[typat]) pspin -> r_max[typat] = pspin -> rrf[typat];
                }
                if (strcmp(str, "h11f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[0][0][3][typat]);
                }
                if (strcmp(str, "h22f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[1][1][3][typat]);
                }
                if (strcmp(str, "h33f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> h[2][2][3][typat]);
                }
                if (strcmp(str, "k11f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k11f[typat]);
                }
                if (strcmp(str, "k22f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k22f[typat]);
                }
                if (strcmp(str, "k33f=") == 0) {
                    fscanf(f2, "%lf", & pspin -> k33f[typat]);
                }
                if (strcmp(str, "COMMENT") == 0) {
                    stop2 = 1;
                }
            }
            printf("PSEUDOPOTENTIAL %d\n", typat + 1);
            printf("   %lf  %lf \n", Z, Z_ion[typat]);
            printf("   rloc= %lf \n", pspin -> rloc[typat]);
            printf("   cc1 = %lf; cc2 = %lf; cc3 = %lf; cc4 = %lf \n", pspin -> cc1[typat], pspin -> cc2[typat], pspin -> cc3[typat], pspin -> cc4[typat]);
            printf("   rrs = %lf; h11s= %lf; h22s= %lf; h33s= %lf \n", pspin -> rrs[typat], pspin -> h[0][0][0][typat], pspin -> h[1][1][0][typat], pspin -> h[2][2][0][typat]);
            printf("   rrp = %lf; h11p= %lf; h22p= %lf; h33p= %lf \n", pspin -> rrp[typat], pspin -> h[0][0][1][typat], pspin -> h[1][1][1][typat], pspin -> h[2][2][1][typat]);
            printf("   rrp = %lf; k11p= %lf; k22p= %lf; k33p= %lf \n", pspin -> rrp[typat], pspin -> k11p[typat], pspin -> k22p[typat], pspin -> k33p[typat]);
            printf("   rrd = %lf; h11d= %lf; h22d= %lf; h33d= %lf \n", pspin -> rrd[typat], pspin -> h[0][0][2][typat], pspin -> h[1][1][2][typat], pspin -> h[2][2][2][typat]);
            printf("   rrd = %lf; k11d= %lf; k22d= %lf; k33d= %lf \n", pspin -> rrd[typat], pspin -> k11d[typat], pspin -> k22d[typat], pspin -> k33d[typat]);
            printf("   rrf = %lf; h11f= %lf; h22f= %lf; h33f= %lf \n", pspin -> rrf[typat], pspin -> h[0][0][3][typat], pspin -> h[1][1][3][typat], pspin -> h[2][2][3][typat]);
            printf("   rrf = %lf; k11f= %lf; k22f= %lf; k33f= %lf \n", pspin -> rrf[typat], pspin -> k11f[typat], pspin -> k22f[typat], pspin -> k33f[typat]);
            printf(" l_max = %d\n", pspin -> l_max[typat]);
            printf(" r_max = %lf\n", pspin -> r_max[typat]);
            printf("\n");
            pspin -> h[0][1][0][typat] = -0.5 * pow(3.0 / 5.0, 0.5) * pspin -> h[1][1][0][typat];
            pspin -> h[1][0][0][typat] = pspin -> h[0][1][0][typat];
            pspin -> h[0][2][0][typat] = 0.5 * pow(5.0 / 21.0, 0.5) * pspin -> h[2][2][0][typat];
            pspin -> h[2][0][0][typat] = pspin -> h[0][2][0][typat];
            pspin -> h[1][2][0][typat] = -0.5 * pow(100.0 / 63.0, 0.5) * pspin -> h[2][2][0][typat];
            pspin -> h[2][1][0][typat] = pspin -> h[1][2][0][typat];
            pspin -> h[0][1][1][typat] = -0.5 * pow(5.0 / 7.0, 0.5) * pspin -> h[1][1][1][typat];
            pspin -> h[1][0][1][typat] = pspin -> h[0][1][1][typat];
            pspin -> h[0][2][1][typat] = (1.0 / 6.0) * pow(35.0 / 11.0, 0.5) * pspin -> h[2][2][1][typat];
            pspin -> h[2][0][1][typat] = pspin -> h[0][2][1][typat];
            pspin -> h[1][2][1][typat] = -(1.0 / 6.0) * 14 * pow(1.0 / 11.0, 0.5) * pspin -> h[2][2][1][typat];
            pspin -> h[2][1][1][typat] = pspin -> h[1][2][1][typat];
            pspin -> h[0][1][2][typat] = -(0.5) * pow(7.0 / 9.0, 0.5) * pspin -> h[1][1][2][typat];
            pspin -> h[1][0][2][typat] = pspin -> h[0][1][2][typat];
            pspin -> h[0][2][2][typat] = (0.5) * pow(63.0 / 143.0, 0.5) * pspin -> h[2][2][2][typat];
            pspin -> h[2][0][2][typat] = pspin -> h[0][2][2][typat];
            pspin -> h[1][2][2][typat] = -(0.5) * 18.0 * pow(1.0 / 143.0, 0.5) * pspin -> h[2][2][2][typat];
            pspin -> h[2][1][2][typat] = pspin -> h[1][2][2][typat];
            typat++;
        }
        if ((strcmp(str, "pspini:") == 0) && (ixc == 11)) {
            finish_line(f2);
            finish_line(f2);
            finish_line(f2);
            fscanf(f2, "%s", str);
            fscanf(f2, "%lf %lf", & Z, & Z_ion[typat]);
            Z_values[typat] = Z;
            pspin -> r_max[typat] = 0.0;
            stop2 = 0;
            l = -1;
            while (stop2 == 0) {
                check = fscanf(f2, "%s", str);
                if (strcmp(str, "angular") == 0) l++;
                if (strcmp(str, "rloc=") == 0) {
                    fscanf(f2, "%lf", & pspin -> rloc[typat]);
                }
                if (strcmp(str, "cc(1:1)=") == 0) {
                    fscanf(f2, " %lf", & pspin -> cc1[typat]);
                }
                if (strcmp(str, "cc2") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc2[typat]);
                }
                if (strcmp(str, "cc3") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc3[typat]);
                }
                if (strcmp(str, "cc4") == 0) {
                    fscanf(f2, " = %lf", & pspin -> cc4[typat]);
                }
                if (strcmp(str, "r(l)") == 0) {
                    if (l == 0) fscanf(f2, " = %lf", & pspin -> rrs[typat]);
                    if (l == 1) fscanf(f2, " = %lf", & pspin -> rrp[typat]);
                    if (l == 2) fscanf(f2, " = %lf", & pspin -> rrd[typat]);
                    if (l == 3) fscanf(f2, " = %lf", & pspin -> rrf[typat]);
                }
                if (strcmp(str, "h11,") == 0) {
                    fscanf(f2, "%s", str);
                    fscanf(f2, "%s", str);
                    if (l == 0) fscanf(f2, " = %lf %lf %lf", & pspin -> h[0][0][0][typat], & pspin -> h[0][1][0][typat], & pspin -> h[0][2][0][typat]);
                    if (l == 1) fscanf(f2, " = %lf %lf %lf", & pspin -> h[0][0][1][typat], & pspin -> h[0][1][1][typat], & pspin -> h[0][2][1][typat]);
                    if (l == 2) fscanf(f2, " = %lf %lf %lf", & pspin -> h[0][0][2][typat], & pspin -> h[0][1][2][typat], & pspin -> h[0][2][2][typat]);
                    if (l == 3) fscanf(f2, " = %lf %lf %lf", & pspin -> h[0][0][3][typat], & pspin -> h[0][1][3][typat], & pspin -> h[0][2][3][typat]);
                }
                if (strcmp(str, "h22,") == 0) {
                    fscanf(f2, "%s", str);
                    if (l == 0) fscanf(f2, " = %lf %lf", & pspin -> h[1][1][0][typat], & pspin -> h[1][2][0][typat]);
                    if (l == 1) fscanf(f2, " = %lf %lf", & pspin -> h[1][1][1][typat], & pspin -> h[1][2][1][typat]);
                    if (l == 2) fscanf(f2, " = %lf %lf", & pspin -> h[1][1][2][typat], & pspin -> h[1][2][2][typat]);
                    if (l == 3) fscanf(f2, " = %lf %lf", & pspin -> h[1][1][3][typat], & pspin -> h[1][2][3][typat]);
                }
                if (strcmp(str, "h33") == 0) {
                    if (l == 0) fscanf(f2, " = %lf", & pspin -> h[2][2][0][typat]);
                    if (l == 1) fscanf(f2, " = %lf", & pspin -> h[2][2][1][typat]);
                    if (l == 2) fscanf(f2, " = %lf", & pspin -> h[2][2][2][typat]);
                    if (l == 3) fscanf(f2, " = %lf", & pspin -> h[2][2][3][typat]);
                }
                if (strcmp(str, "k11,") == 0) {
                    fscanf(f2, "%s", str);
                    fscanf(f2, "%s", str);
                    //                    if(l==0) fscanf(f2,"= %lf %lf %lf",&pspin->k11s[typat],&pspin->k12s[typat],&pspin->k13s[typat]);
                    if (l == 1) fscanf(f2, " = %lf %lf %lf", & pspin -> k11p[typat], & pspin -> k12p[typat], & pspin -> k13p[typat]);
                    if (l == 2) fscanf(f2, " = %lf %lf %lf", & pspin -> k11d[typat], & pspin -> k12d[typat], & pspin -> k13d[typat]);
                    if (l == 3) fscanf(f2, " = %lf %lf %lf", & pspin -> k11f[typat], & pspin -> k12f[typat], & pspin -> k13f[typat]);
                }
                if (strcmp(str, "k22,") == 0) {
                    fscanf(f2, "%s", str);
                    //                    if(l==0) fscanf(f2,"= %lf %lf",&pspin->k22s[typat],&pspin->k23s[typat]);
                    if (l == 1) fscanf(f2, " = %lf %lf", & pspin -> k22p[typat], & pspin -> k23p[typat]);
                    if (l == 2) fscanf(f2, " = %lf %lf", & pspin -> k22d[typat], & pspin -> k23d[typat]);
                    if (l == 3) fscanf(f2, " = %lf %lf", & pspin -> k22f[typat], & pspin -> k23f[typat]);
                }
                if (strcmp(str, "k33") == 0) {
                    fscanf(f2, "%s", str);
                    //                   if(l==0) fscanf(f2,"= %lf",&pspin->k33s[typat]);
                    if (l == 1) fscanf(f2, " = %lf", & pspin -> k33p[typat]);
                    if (l == 2) fscanf(f2, " = %lf", & pspin -> k33d[typat]);
                    if (l == 3) fscanf(f2, " = %lf", & pspin -> k33f[typat]);
                }
                if (strcmp(str, "COMMENT") == 0) {
                    stop2 = 1;
                }
            }
            if (pspin -> rrs[typat] > 0) pspin -> l_max[typat] = 1;
            if (pspin -> rrp[typat] > 0) pspin -> l_max[typat] = 2;
            if (pspin -> rrd[typat] > 0) pspin -> l_max[typat] = 3;
            if (pspin -> rrf[typat] > 0) pspin -> l_max[typat] = 4;
            if (pspin -> rrs[typat] > pspin -> r_max[typat]) pspin -> r_max[typat] = pspin -> rrs[typat];
            if (pspin -> rrp[typat] > pspin -> r_max[typat]) pspin -> r_max[typat] = pspin -> rrp[typat];
            if (pspin -> rrd[typat] > pspin -> r_max[typat]) pspin -> r_max[typat] = pspin -> rrd[typat];
            if (pspin -> rrf[typat] > pspin -> r_max[typat]) pspin -> r_max[typat] = pspin -> rrf[typat];
            printf("PSEUDOPOTENTIAL %d\n", typat + 1);
            printf("   %lf  %lf \n", Z, Z_ion[typat]);
            printf("   rloc= %lf \n", pspin -> rloc[typat]);
            printf("   cc1 = %lf; cc2 = %lf; cc3 = %lf; cc4 = %lf \n", pspin -> cc1[typat], pspin -> cc2[typat], pspin -> cc3[typat], pspin -> cc4[typat]);
            printf("   rrs = %lf; h11s= %lf; h22s= %lf; h33s= %lf \n", pspin -> rrs[typat], pspin -> h[0][0][0][typat], pspin -> h[1][1][0][typat], pspin -> h[2][2][0][typat]);
            printf("   rrp = %lf; h11p= %lf; h22p= %lf; h33p= %lf \n", pspin -> rrp[typat], pspin -> h[0][0][1][typat], pspin -> h[1][1][1][typat], pspin -> h[2][2][1][typat]);
            printf("   rrp = %lf; k11p= %lf; k22p= %lf; k33p= %lf \n", pspin -> rrp[typat], pspin -> k11p[typat], pspin -> k22p[typat], pspin -> k33p[typat]);
            printf("   rrd = %lf; h11d= %lf; h22d= %lf; h33d= %lf \n", pspin -> rrd[typat], pspin -> h[0][0][2][typat], pspin -> h[1][1][2][typat], pspin -> h[2][2][2][typat]);
            printf("   rrd = %lf; k11d= %lf; k22d= %lf; k33d= %lf \n", pspin -> rrd[typat], pspin -> k11d[typat], pspin -> k22d[typat], pspin -> k33d[typat]);
            printf("   rrf = %lf; h11f= %lf; h22f= %lf; h33f= %lf \n", pspin -> rrf[typat], pspin -> h[0][0][3][typat], pspin -> h[1][1][3][typat], pspin -> h[2][2][3][typat]);
            printf("   rrf = %lf; k11f= %lf; k22f= %lf; k33f= %lf \n", pspin -> rrf[typat], pspin -> k11f[typat], pspin -> k22f[typat], pspin -> k33f[typat]);
            printf(" l_max = %d\n", pspin -> l_max[typat]);
            printf(" r_max = %lf\n", pspin -> r_max[typat]);
            printf("\n");
            //               pspin->h[0][1][0][typat]=-0.5*pow(3.0/5.0,0.5)*pspin->h[1][1][0][typat];
            pspin -> h[1][0][0][typat] = pspin -> h[0][1][0][typat];
            //               pspin->h[0][2][0][typat]=0.5*pow(5.0/21.0,0.5)*pspin->h[2][2][0][typat];
            pspin -> h[2][0][0][typat] = pspin -> h[0][2][0][typat];
            //               pspin->h[1][2][0][typat]=-0.5*pow(100.0/63.0,0.5)*pspin->h[2][2][0][typat];
            pspin -> h[2][1][0][typat] = pspin -> h[1][2][0][typat];
            //               pspin->h[0][1][1][typat]=-0.5*pow(5.0/7.0,0.5)*pspin->h[1][1][1][typat];
            pspin -> h[1][0][1][typat] = pspin -> h[0][1][1][typat];
            //              pspin->h[0][2][1][typat]=(1.0/6.0)*pow(35.0/11.0,0.5)*pspin->h[2][2][1][typat];
            pspin -> h[2][0][1][typat] = pspin -> h[0][2][1][typat];
            //               pspin->h[1][2][1][typat]=-(1.0/6.0)*14*pow(1.0/11.0,0.5)*pspin->h[2][2][1][typat];
            pspin -> h[2][1][1][typat] = pspin -> h[1][2][1][typat];
            //               pspin->h[0][1][2][typat]=-(0.5)*pow(7.0/9.0,0.5)*pspin->h[1][1][2][typat];
            pspin -> h[1][0][2][typat] = pspin -> h[0][1][2][typat];
            //               pspin->h[0][2][2][typat]=(0.5)*pow(63.0/143.0,0.5)*pspin->h[2][2][2][typat];
            pspin -> h[2][0][2][typat] = pspin -> h[0][2][2][typat];
            //               pspin->h[1][2][2][typat]=-(0.5)*18.0*pow(1.0/143.0,0.5)*pspin->h[2][2][2][typat];
            pspin -> h[2][1][2][typat] = pspin -> h[1][2][2][typat];
            printf(" h11s h12s h13s = %lf %lf %lf \n", pspin -> h[0][0][0][typat], pspin -> h[0][1][0][typat], pspin -> h[0][2][0][typat]);
            printf("      h22s h23s =     %lf %lf \n", pspin -> h[1][1][0][typat], pspin -> h[1][2][0][typat]);
            printf("           h33s =     %lf %lf \n", pspin -> h[2][2][0][typat]);
            printf(" h11p h12p h13p = %lf %lf %lf \n", pspin -> h[0][0][1][typat], pspin -> h[0][1][1][typat], pspin -> h[0][2][1][typat]);
            printf("      h22p h23p =     %lf %lf \n", pspin -> h[1][1][1][typat], pspin -> h[1][2][1][typat]);
            printf("           h33p =     %lf %lf \n", pspin -> h[2][2][1][typat]);
            printf(" h11d h12d h13d = %lf %lf %lf \n", pspin -> h[0][0][2][typat], pspin -> h[0][1][2][typat], pspin -> h[0][2][2][typat]);
            printf("      h22d h23d =     %lf %lf \n", pspin -> h[1][1][2][typat], pspin -> h[1][2][2][typat]);
            printf("           h33d =     %lf %lf \n", pspin -> h[2][2][2][typat]);
            printf(" h11f h12f h13f = %lf %lf %lf \n", pspin -> h[0][0][3][typat], pspin -> h[0][1][3][typat], pspin -> h[0][2][3][typat]);
            printf("      h22f h23f =     %lf %lf \n", pspin -> h[1][1][3][typat], pspin -> h[1][2][3][typat]);
            printf("           h33f =     %lf %lf \n", pspin -> h[2][2][3][typat]);
            typat++;
        }
    }
    for (j = 0; j < typat; j++) {
        if (l_max < pspin -> l_max[j]) l_max = pspin -> l_max[j];
    }
    printf("  Overall l_max = %d", l_max);
    fclose(f2);
}

double projp(int l, int i, double g1, struct psp_file * pspin, int typat, double cellvolume) {
    //                                                  p1=projp(l,i,g1,pspin,atom_type,cellvolume);
    double p = 0.0;
    typat--;
    if ((l == 0) && (i == 0)) {
        p = 4.0 * pow(2.0 * pow(pspin -> rrs[typat], 3.0), 0.5) * PI_5_4 / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrs[typat], 2.0)));
        //	       printf("p = %lf   rrs = %lf g1 = %lf atom_type = %d,  cell_volume = %lf\n",p, pspin->rrs[0],g1,typat,cellvolume);
        //               exit(0);
    }
    if ((l == 0) && (i == 1)) {
        p = 8.0 * pow(2.0 * pow(pspin -> rrs[typat], 3.0) / 15.0, 0.5) * PI_5_4 * (3.0 - pow(g1 * pspin -> rrs[typat], 2.0)) / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrs[typat], 2.0)));
    }
    if ((l == 0) && (i == 2)) {
        p = 16.0 * pow(2.0 * pow(pspin -> rrs[typat], 3.0) / 105.0, 0.5) * PI_5_4 * (15.0 - 10.0 * pow(g1 * pspin -> rrs[typat], 2.0) + pow(g1 * pspin -> rrs[typat], 4.0)) / (3.0 * pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrs[typat], 2.0)));
    }
    if ((l == 1) && (i == 0)) {
        p = 8.0 * pow(pow(pspin -> rrp[typat], 5.0) / 3.0, 0.5) * PI_5_4 * (g1) / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrp[typat], 2.0)));
    }
    if ((l == 1) && (i == 1)) {
        p = 16.0 * pow(pow(pspin -> rrp[typat], 5.0) / 105.0, 0.5) * PI_5_4 * (g1) * (5.0 - pow(g1 * pspin -> rrp[typat], 2.0)) / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrp[typat], 2.0)));
    }
    if ((l == 1) && (i == 2)) {
        p = 32.0 * pow(pow(pspin -> rrp[typat], 5.0) / 1155.0, 0.5) * PI_5_4 * (g1) * (35.0 - 14.0 * pow(g1 * pspin -> rrp[typat], 2.0) + pow(g1 * pspin -> rrp[typat], 4.0)) / (3.0 * pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrp[typat], 2.0)));
    }
    if ((l == 2) && (i == 0)) {
        p = 8.0 * pow(2 * pow(pspin -> rrd[typat], 7.0) / 15.0, 0.5) * PI_5_4 * (g1 * g1) / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrd[typat], 2.0)));
    }
    if ((l == 2) && (i == 1)) {
        p = 16.0 * pow(2 * pow(pspin -> rrd[typat], 7.0) / 105.0, 0.5) * PI_5_4 * (g1 * g1) * (7.0 - pow(g1 * pspin -> rrd[typat], 2.0)) / (3.0 * pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrd[typat], 2.0)));
    }
    if ((l == 3) && (i == 0)) {
        p = 16.0 * pow(pow(pspin -> rrf[typat], 9.0) / 105.0, 0.5) * PI_5_4 * (g1 * g1 * g1) / (pow(cellvolume, 0.5) * exp(0.5 * pow(g1 * pspin -> rrf[typat], 2.0)));
    }
    //        if(p<0.0) printf("YIKES!");
    return (p);
}

void read_calc_nonlocal(char filename[200], char outfile[200], struct wfk_file * wfkin, struct psp_file * pspin, double occ_min, int dtset) {
    FILE * f2;
    FILE * f3;
    int k;
    char codvsn[110];
    char title[132];
    char test;
    int headform;
    int fform;
    int bandtot;
    int date;
    int intxc;
    int ixc;
    int natom;
    int atomno, pw1, pw2;
    int ngfftx;
    int ngffty;
    int ngfftz;
    int nkpt;
    int npsp;
    int j = 0;
    int i, m, l, ka, kb, kc;
    int stop = 0;
    int nspden;
    int nsppol;
    int ntypat;
    int occopt;
    int pertcase;
    int usepaw;
    double ecut;
    double ecutdg;
    double ecutsm;
    double ecut_eff;
    double qptnx;
    double qptny;
    double qptnz;
    double rprimd_ax;
    double rprimd_ay;
    double rprimd_az;
    double rprimd_bx;
    double rprimd_by;
    double rprimd_bz;
    double rprimd_cx;
    double rprimd_cy;
    double rprimd_cz;
    double ax_star;
    double ay_star;
    double az_star;
    double bx_star;
    double by_star;
    double bz_star;
    double cx_star;
    double cy_star;
    double cz_star;
    double ox, oy, oz;
    double theta1, phi1, g1, theta3;
    double theta2, phi2, g2;
    double stmbias;
    double tphysel;
    double tsmear;
    double znuclpsp;
    double zionpsp;
    double costheta1;
    double cos_theta2;
    double cos_table[10000000];
    double sin_table[10000000];
    double theta_step = 0.000001;
    int pspso;
    int pspdat;
    int pspcod;
    int pspxc;
    int lmn_size;
    int atom_type, atom;
    int usewvl;
    int istwfk[KPOINTS_MAX];
    int nband[BANDS_MAX];
    int npwarr[KPOINTS_MAX];
    int so_psp[NATOMS_MAX];
    int symafm[MAX_SYM];
    double kpt[3][KPOINTS_MAX];
    double occ[TBANDS_MAX];
    double cgt[2][WAVES_MAX][2];
    double cg_band[WAVES_MAX][BANDS_MAX][2];
    double cg_beta[WAVES_MAX][BANDS_MAX][2];
    double Int_PjYlm_RE[4][8][4][NATOMS_MAX]; /* l,m,j,atom */
    double Int_PjYlm_IM[4][8][4][NATOMS_MAX]; /* l,m,j,atom */
    double znucltypat[NATOMS_MAX];
    double wtk[KPOINTS_MAX];
    double residm;
    double x, y, z;
    double etotal, fermie;
    double Enonlocal_temp;
    /* p1=projp(l,i,g[pw1],pspin,atom_type,cellvolume); */
    double ga1, gb1, gc1;
    double gx1, gy1, gz1;
    double ga2, gb2, gc2;
    double gx2, gy2, gz2;
    double xred[3][NATOMS_MAX];
    double php, yphpy;
    double yphpy_l[4];
    int pw_counter;
    double frac_pw;
    int npw;
    int pw;
    int nspinor;
    int spin;
    int nband_temp;
    int noccband;
    int kptno;
    int kx, ky, kz, jx, jy, jz, k_ngx, k_ngy, k_ngz;
    int offsetRE, offsetIM;
    int band;
    int i_on[NTYPES_MAX][5][5];
    int l_on[NTYPES_MAX][5];
    double Dga, Dgb, Dgc;
    double normalization;
    double eigen, occ_temp, cg;
    double cellvolume, p1, p2;
    gsl_complex c1;
    gsl_complex c2;
    gsl_complex c1_star;
    gsl_complex c1c2;
    gsl_complex atom_phase;
    gsl_complex gkk;
    double atom_phase0;
    double atom_phase_rad;
    double Dphi, cos_mdphi_times_2[4];
    gsl_complex vg1g2;
    double vg1g2_real, vg1g2_real_temp;
    double vg1g2_real_l[NTYPES_MAX][4];
    double c1c2atomphase_real;
    double c1c2_RE[BANDS_MAX], c1c2_IM[BANDS_MAX], gkk_real;
    double atomphase_RE, atomphase_IM;
    gsl_complex vg1g2_lm;
    gsl_complex vg1g2_temp;
    int l_max_type[NTYPES_MAX];
    int i_max_type[NTYPES_MAX][4];
    int table_entry, table_entry1, table_entry2, table_entry3;
    double ph[4][4][5][NTYPES_MAX];
    fftw_complex * psi_R, * psi_G;
    double * psi_R_ptr = (double * ) psi_R;
    double * psi_G_ptr = (double * ) psi_G;
    fftw_complex * gradx_R, * gradx_G;
    fftw_complex * grady_R, * grady_G;
    fftw_complex * gradz_R, * gradz_G;
    double * gradx_R_ptr = (double * ) gradx_R;
    double * grady_R_ptr = (double * ) grady_R;
    double * gradz_R_ptr = (double * ) gradz_R;
    double * gradx_G_ptr = (double * ) gradx_G;
    double * grady_G_ptr = (double * ) grady_G;
    double * gradz_G_ptr = (double * ) gradz_G;
    fftw_plan fftplan_psi;
    fftw_plan fftplan_psi_inv;
    fftw_plan fftplan_psi_inv2;
    fftw_plan fftplan_gradx;
    fftw_plan fftplan_grady;
    fftw_plan fftplan_gradz;
    double xc, yc, zc, voxcenter_x, voxcenter_y, voxcenter_z, theta, phi;
    double sinmphi[4], cosmphi[4], dist;
    double xf, yf, zf;
    /* END_DECLARATIONS */
    printf("Starting non-local calculation\n");
    theta1 = 0;
    for (j = 0; j < 9999999; j++) {
        cos_table[j] = cos(TWO_PI * theta1);
        sin_table[j] = sin(TWO_PI * theta1);
        theta1 += theta_step;
    }
    theta1 = 0;
    for (j = 0; j < NATOMS_MAX; j++) {
        for (l = 0; l < 5; l++) {
            nonlocalE_byatom[dtset][j][l] = 0.0;
        }
    }
    f2 = fopen(filename, "rb+");
    if (f2 == NULL) {
        printf("%s not found.\n", filename);
        exit(0);
    }
    /*   READ HEADER OF WFK FILE    */
    j = 0;
    fread( & j, sizeof(int), 1, f2);
    j = fread(codvsn, sizeof(char), 6, f2);
    fread( & headform, sizeof(int), 1, f2);
    fread( & fform, sizeof(int), 1, f2);
    printf("%s %d %d\n", codvsn, headform, fform);
    fread( & j, sizeof(int), 1, f2);
    fread( & j, sizeof(int), 1, f2);
    fread( & bandtot, sizeof(int), 1, f2);
    /*    int date;*/
    fread( & date, sizeof(int), 1, f2);
    printf("%d %d\n", bandtot, date);
    /*    int intxc;  */
    fread( & intxc, sizeof(int), 1, f2);
    /*    int ixc;*/
    fread( & ixc, sizeof(int), 1, f2);
    /*    int natom;*/
    fread( & natom, sizeof(int), 1, f2);
    if (natom > NATOMS_MAX) {
        printf("natom > NATOMS_MAX!\n");
        exit(0);
    }
    /*    int ngfftx;*/
    fread( & ngfftx, sizeof(int), 1, f2);
    /*    int ngffty;*/
    fread( & ngffty, sizeof(int), 1, f2);
    /*    int ngfftz;*/
    fread( & ngfftz, sizeof(int), 1, f2);
    printf("ngfft = %d x %d x %d\n", ngfftx, ngffty, ngfftz);
    /*    int nkpt;*/
    den1.NGX = ngfftx + 1;
    den1.NGY = ngffty + 1;
    den1.NGZ = ngfftz + 1;
    printf("Allocating space of %d x %d x %d grids in real and k space.\n", ngfftx, ngffty, ngfftz);
    psi_R = (fftw_complex * ) fftw_malloc(sizeof(fftw_complex) * ngfftx * ngffty * ngfftz);
    psi_G = (fftw_complex * ) fftw_malloc(sizeof(fftw_complex) * ngfftx * ngffty * ngfftz);
    fftplan_psi = fftw_plan_dft_3d(ngfftx, ngffty, ngfftz, psi_G, psi_R, FFTW_FORWARD, FFTW_ESTIMATE);
    fftplan_psi_inv = fftw_plan_dft_3d(ngfftx, ngffty, ngfftz, psi_R, psi_G, FFTW_FORWARD, FFTW_ESTIMATE);
    fftplan_psi_inv2 = fftw_plan_dft_3d(ngfftx, ngffty, ngfftz, psi_G, psi_R, FFTW_BACKWARD, FFTW_ESTIMATE);
    psi_R_ptr = (double * ) psi_R;
    psi_G_ptr = (double * ) psi_G;
    fread( & nkpt, sizeof(int), 1, f2);
    if (nkpt > KPOINTS_MAX) {
        printf("nkpt > KPOINTS_MAX!\n");
        exit(0);
    }
    /*    int nspden;*/
    fread( & nspden, sizeof(int), 1, f2);
    /*    int nspinor;*/
    fread( & nspinor, sizeof(int), 1, f2);
    /*    int nsppol;*/
    if (nspinor > 1) printf("  Detected nspinor = %d\n", nspinor);
    fread( & nsppol, sizeof(int), 1, f2);
    /*    int nsym;*/
    fread( & nsym, sizeof(int), 1, f2);
    /*    int npsp;*/
    fread( & npsp, sizeof(int), 1, f2);
    /*    int ntypat;*/
    fread( & ntypat, sizeof(int), 1, f2);
    /*    int occopt;*/
    fread( & occopt, sizeof(int), 1, f2);
    /*    int pertcase;*/
    fread( & pertcase, sizeof(int), 1, f2);
    /*    int usepaw;*/
    fread( & usepaw, sizeof(int), 1, f2);
    /*    double ecut;*/
    fread( & ecut, sizeof(double), 1, f2);
    /*     double ecutdg;*/
    fread( & ecutdg, sizeof(double), 1, f2);
    /*    double ecutsm;*/
    fread( & ecutsm, sizeof(double), 1, f2);
    /*    double ecut_eff;*/
    fread( & ecut_eff, sizeof(double), 1, f2);
    /*    double qptnx;*/
    fread( & qptnx, sizeof(double), 1, f2);
    /*    double qptny;*/
    fread( & qptny, sizeof(double), 1, f2);
    /*    double qptnz;*/
    fread( & qptnz, sizeof(double), 1, f2);
    /*    double rprimd_ax;*/
    fread( & rprimd_ax, sizeof(double), 1, f2);
    /*    double rprimd_ay;*/
    fread( & rprimd_ay, sizeof(double), 1, f2);
    /*    double rprimd_az;*/
    fread( & rprimd_az, sizeof(double), 1, f2);
    /*    double rprimd_bx;*/
    fread( & rprimd_bx, sizeof(double), 1, f2);
    /*    double rprimd_by;*/
    fread( & rprimd_by, sizeof(double), 1, f2);
    /*    double rprimd_bz;*/
    fread( & rprimd_bz, sizeof(double), 1, f2);
    /*    double rprimd_cx;*/
    fread( & rprimd_cx, sizeof(double), 1, f2);
    /*    double rprimd_cy;*/
    fread( & rprimd_cy, sizeof(double), 1, f2);
    /*    double rprimd_cz;*/
    fread( & rprimd_cz, sizeof(double), 1, f2);
    /*    double stmbias;*/
    den0.cella_x = rprimd_ax;
    den0.cella_y = rprimd_ay;
    den0.cella_z = rprimd_az;
    den0.cellb_x = rprimd_bx;
    den0.cellb_y = rprimd_by;
    den0.cellb_z = rprimd_bz;
    den0.cellc_x = rprimd_cx;
    den0.cellc_y = rprimd_cy;
    den0.cellc_z = rprimd_cz;
    cellvolume = (rprimd_ax * (rprimd_by * rprimd_cz - rprimd_bz * rprimd_cy) - rprimd_ay * (rprimd_bx * rprimd_cz - rprimd_bz * rprimd_cx) + rprimd_az * (rprimd_bx * rprimd_cy - rprimd_by * rprimd_cx));
    printf(" cella = %lf %lf %lf \n", rprimd_ax, rprimd_ay, rprimd_az);
    printf(" cellb = %lf %lf %lf \n", rprimd_bx, rprimd_by, rprimd_bz);
    printf(" cellc = %lf %lf %lf \n\n", rprimd_cx, rprimd_cy, rprimd_cz);
    printf(" cell volume = %lf\n\n", cellvolume);
    ax_star = 2 * PI * (rprimd_by * rprimd_cz - rprimd_cy * rprimd_bz) / cellvolume;
    ay_star = -2 * PI * (rprimd_bx * rprimd_cz - rprimd_cx * rprimd_bz) / cellvolume;
    az_star = 2 * PI * (rprimd_bx * rprimd_cy - rprimd_cx * rprimd_by) / cellvolume;
    bx_star = 2 * PI * (rprimd_cy * rprimd_az - rprimd_ay * rprimd_cz) / cellvolume;
    by_star = -2 * PI * (rprimd_cx * rprimd_az - rprimd_ax * rprimd_cz) / cellvolume;
    bz_star = 2 * PI * (rprimd_cx * rprimd_ay - rprimd_ax * rprimd_cy) / cellvolume;
    cx_star = 2 * PI * (rprimd_ay * rprimd_bz - rprimd_by * rprimd_az) / cellvolume;
    cy_star = -2 * PI * (rprimd_ax * rprimd_bz - rprimd_bx * rprimd_az) / cellvolume;
    cz_star = 2 * PI * (rprimd_ax * rprimd_by - rprimd_bx * rprimd_ay) / cellvolume;
    printf(" cella* = %lf %lf %lf \n", ax_star, ay_star, az_star);
    printf(" cellb* = %lf %lf %lf \n", bx_star, by_star, bz_star);
    printf(" cellc* = %lf %lf %lf \n\n", cx_star, cy_star, cz_star);
    printf(" a.a*=%lf, a.b*=%lf, a.c*=%lf\n", rprimd_ax * ax_star + rprimd_ay * ay_star + rprimd_az * az_star, rprimd_ax * bx_star + rprimd_ay * by_star + rprimd_az * bz_star, rprimd_ax * cx_star + rprimd_ay * cy_star + rprimd_az * cz_star);
    printf(" b.a*=%lf, b.b*=%lf, b.c*=%lf\n", rprimd_bx * ax_star + rprimd_by * ay_star + rprimd_bz * az_star, rprimd_bx * bx_star + rprimd_by * by_star + rprimd_bz * bz_star, rprimd_bx * cx_star + rprimd_by * cy_star + rprimd_bz * cz_star);
    printf(" c.a*=%lf, c.b*=%lf, c.c*=%lf\n\n", rprimd_cx * ax_star + rprimd_cy * ay_star + rprimd_cz * az_star, rprimd_cx * bx_star + rprimd_cy * by_star + rprimd_cz * bz_star, rprimd_cx * cx_star + rprimd_cy * cy_star + rprimd_cz * cz_star);

    fread( & stmbias, sizeof(double), 1, f2);
    /*    double tphysel;*/
    fread( & tphysel, sizeof(double), 1, f2);
    /*    double tsmear;*/
    fread( & tsmear, sizeof(double), 1, f2);
    /*    int usewvl;*/
    fread( & usewvl, sizeof(int), 1, f2);
    //    printf("natoms = %d   ecut = %lf  tsmear = %lf  occopt = %d \n",natom,ecut,tsmear, occopt); 

    fread( & j, sizeof(int), 1, f2);
    fread( & j, sizeof(int), 1, f2);
    /*    int istwfk [KPOINTS_MAX]; */
    //    printf("nkpt = %d    nsppol = %d   npsp = %d  ntypat = %d \n",nkpt,nsppol,npsp,ntypat);
    for (j = 0; j < nkpt; j++) {
        fread( & istwfk[j], sizeof(int), 1, f2);
        if (istwfk[j] > 1) {
            printf("\n\nUH OH.  \n\nistwfk value for kpoint %d is %d.  Please manually set to all istwfk values to 1,\n", j, istwfk[j]);
            printf("in _in file and rerun abinit.\n\n", j, istwfk[j]);
            exit(0);
        }
    }
    /*    int nband [BANDS_MAX]; */
    for (j = 0; j < (nkpt * nsppol); j++) {
        fread( & nband[j], sizeof(int), 1, f2);
        if (nband[j] > BANDS_MAX) {
            printf("nband at kpt %d  (%d) > BANDS_MAX.\n", j + 1, nband[j]);
            exit(0);
        }
    }
    /*    int npwarr [KPOINTS_MAX]; */
    for (j = 0; j < (nkpt); j++) {
        fread( & npwarr[j], sizeof(int), 1, f2);
    }
    /*    int so_psp [NATOMS_MAX]; */
    printf("\n");
    for (j = 0; j < (npsp); j++) {
        fread( & so_psp[j], sizeof(int), 1, f2);
        /* printf("PSEUDOPOTENTIAL %d\n",j+1);
	    printf("   rloc= %lf \n",pspin->rloc[j]);
	    printf("   cc1 = %lf; cc2 = %lf; cc3 = %lf; cc4 = %lf \n",pspin->cc1[j],pspin->cc2[j],pspin->cc3[j],pspin->cc4[j]);
	    printf("   rrs = %lf; h11s= %lf; h22s= %lf; h33s= %lf \n",pspin->rrs[j],pspin->h[0][0][0][j],pspin->h[1][1][0][j],pspin->h[2][2][0][j]);
	    printf("   rrp = %lf; h11p= %lf; h22p= %lf; h33p= %lf \n",pspin->rrp[j],pspin->h[0][0][1][j],pspin->h[1][1][1][j],pspin->h[2][2][1][j]);
	    printf("   rrp = %lf; k11p= %lf; k22p= %lf; k33p= %lf \n",pspin->rrp[j],pspin->k11p[j],pspin->k22p[j],pspin->k33p[j]);
	    printf("   rrd = %lf; h11d= %lf; h22d= %lf; h33d= %lf \n",pspin->rrd[j],pspin->h[0][0][2][j],pspin->h[1][1][2][j],pspin->h[2][2][2][j]);
	    printf("   rrd = %lf; k11d= %lf; k22d= %lf; k33d= %lf \n",pspin->rrd[j],pspin->k11d[j],pspin->k22d[j],pspin->k33d[j]);
	    printf("\n");  */
    }
    /*    int symafm [MAX_SYM]; */
    for (j = 0; j < (nsym); j++) {
        fread( & symafm[j], sizeof(int), 1, f2);
    }
    /*    int symrel [3][3][MAX_SYM]; */
    for (j = 0; j < (nsym); j++) {
        fread( & symrel[0][0][j], sizeof(int), 1, f2);
        fread( & symrel[1][0][j], sizeof(int), 1, f2);
        fread( & symrel[2][0][j], sizeof(int), 1, f2);
        fread( & symrel[0][1][j], sizeof(int), 1, f2);
        fread( & symrel[1][1][j], sizeof(int), 1, f2);
        fread( & symrel[2][1][j], sizeof(int), 1, f2);
        fread( & symrel[0][2][j], sizeof(int), 1, f2);
        fread( & symrel[1][2][j], sizeof(int), 1, f2);
        fread( & symrel[2][2][j], sizeof(int), 1, f2);
    }
    /*    int typat [NATOMS_MAX]; */
    for (j = 0; j < (natom); j++) {
        fread( & typat[j], sizeof(int), 1, f2);
        //	 printf("%d\n",typat[j]);         
    }
    /*    double kpt [3][KPOINTS_MAX]; */
    for (j = 0; j < (nkpt); j++) {
        fread( & kpt[0][j], sizeof(double), 1, f2);
        fread( & kpt[1][j], sizeof(double), 1, f2);
        fread( & kpt[2][j], sizeof(double), 1, f2);
        //         printf("kpoint %d:  %lf %lf %lf \n",j+1,kpt[0][j],kpt[1][j],kpt[2][j]);
    }
    /*    double occ(TBANDS_MAX); */
    for (j = 0; j < (bandtot); j++) {
        fread( & occ[0], sizeof(double), 1, f2);
    }
    /*    double tnons [3][MAX_SYM]; */
    for (j = 0; j < (nsym); j++) {
        fread( & tnons[0][j], sizeof(double), 1, f2);
        fread( & tnons[1][j], sizeof(double), 1, f2);
        fread( & tnons[2][j], sizeof(double), 1, f2);
    }
    /*    double znucltypat[NATOMS_MAX];  */
    for (j = 0; j < (ntypat); j++) {
        fread( & znucltypat[j], sizeof(double), 1, f2);
    }
    /*    double wtk[KPOINTS_MAX]; */
    for (j = 0; j < (nkpt); j++) {
        fread( & wtk[j], sizeof(double), 1, f2);
    }
    fread( & j, sizeof(int), 1, f2);
    for (k = 0; k < (npsp); k++) {
        fread( & j, sizeof(int), 1, f2);
        fread(title, sizeof(char), 132, f2);
        //         printf("     %s\n",title);
        fread( & znuclpsp, sizeof(double), 1, f2);
        //          printf("     %lf ",znuclpsp);
        fread( & zionpsp, sizeof(double), 1, f2);
        //          printf("%lf ",zionpsp);
        fread( & pspso, sizeof(int), 1, f2);
        fread( & pspdat, sizeof(int), 1, f2);
        fread( & pspcod, sizeof(int), 1, f2);
        fread( & pspxc, sizeof(int), 1, f2);
        fread( & lmn_size, sizeof(int), 1, f2);
        //          printf("%d %d %d %d %d \n",pspso,pspdat,pspcod,pspxc,lmn_size);
        fread( & j, sizeof(int), 1, f2);
    }
    if (usepaw == 0) {
        fread( & j, sizeof(int), 1, f2);
        fread( & residm, sizeof(double), 1, f2);
        //        printf("     residm = %lf\n",residm);
        for (k = 0; k < natom; k++) {
            fread( & xred[0][k], sizeof(double), 1, f2);
            fread( & xred[1][k], sizeof(double), 1, f2);
            fread( & xred[2][k], sizeof(double), 1, f2);
            printf("     Atom %d:  %lf %lf %lf \n", k + 1, xred[0][k], xred[1][k], xred[2][k]);
            x = xred[0][k];
            y = xred[1][k];
            z = xred[2][k];
            den1.Xcart[k] = (x) * den1.cella_x + (y) * den1.cellb_x + (z) * den1.cellc_x;
            den1.Ycart[k] = (x) * den1.cella_y + (y) * den1.cellb_y + (z) * den1.cellc_y;
            den1.Zcart[k] = (x) * den1.cella_z + (y) * den1.cellb_z + (z) * den1.cellc_z;
        }
        fread( & etotal, sizeof(double), 1, f2);
        fread( & fermie, sizeof(double), 1, f2);
        printf("     Etotal = %lf    FermiE = %lf\n\n", etotal, fermie);
        fread( & j, sizeof(int), 1, f2);
    } else {
        printf("Yikes!  usepaw!=0.  We haven't written code for this case yet. \n");
    }
    /*   HEADER FINISHED - READ WFK DATA, CALCULATE NONLOCAL ENERGY CONTRIBUTIONS   */
    cos_mdphi_times_2[0] = 0.0;
    for (jx = 0; jx < den1.NGX - 1; jx++) {
        for (jy = 0; jy < den1.NGY - 1; jy++) {
            for (jz = 0; jz < den1.NGZ - 1; jz++) {
                den1.grid[jx][jy][jz] = 0.0;
            }
        }
    }
    for (k = 0; k < nsppol; k++) {
        for (kptno = 0; kptno < nkpt; kptno++) {
            printf("Calculating for kpoint %d:  ", kptno + 1);
            fread( & j, sizeof(int), 1, f2);
            fread( & npw, sizeof(int), 1, f2);
            if (npw > WAVES_MAX) {
                printf("npw > WAVES_MAX at k-point %d!\n", j + 1);
                exit(0);
            }
            fread( & nspinor, sizeof(int), 1, f2);
            fread( & nband_temp, sizeof(int), 1, f2);
            fread( & j, sizeof(int), 1, f2);
            fread( & j, sizeof(int), 1, f2);
            for (pw = 0; pw < npw; pw++) {
                fread( & kx, sizeof(int), 1, f2);
                fread( & ky, sizeof(int), 1, f2);
                fread( & kz, sizeof(int), 1, f2);
                wfkin -> kg[pw][0] = kx;
                wfkin -> kg[pw][1] = ky;
                wfkin -> kg[pw][2] = kz;
            }
            printf(" nbands = %d,  nwaves = %d, weight = %lf\n", nband_temp, npw, wtk[kptno]);
            fread( & j, sizeof(int), 1, f2);
            fread( & j, sizeof(int), 1, f2);
            for (band = 0; band < nband_temp; band++) {
                fread( & eigen, sizeof(double), 1, f2);
                wfkin -> eigen[band] = eigen;
            }
            noccband = 0;
            for (band = 0; band < nband_temp; band++) {
                fread( & occ_temp, sizeof(double), 1, f2);
                occ[band] = wtk[kptno] * occ_temp;
                if (occ_temp > occ_min) noccband++;
            }
            fread( & j, sizeof(int), 1, f2);
            for (pw1 = 0; pw1 < npw; pw1++) {
                ga[pw1] = kpt[0][kptno] + wfkin -> kg[pw1][0] * 1.0;
                gb[pw1] = kpt[1][kptno] + wfkin -> kg[pw1][1] * 1.0;
                gc[pw1] = kpt[2][kptno] + wfkin -> kg[pw1][2] * 1.0;
                gx[pw1] = ga[pw1] * ax_star + gb[pw1] * bx_star + gc[pw1] * cx_star;
                gy[pw1] = ga[pw1] * ay_star + gb[pw1] * by_star + gc[pw1] * cy_star;
                gz[pw1] = ga[pw1] * az_star + gb[pw1] * bz_star + gc[pw1] * cz_star;

            }
            for (band = 0; band < nband_temp; band++) {
                fread( & j, sizeof(int), 1, f2);
                //              printf("   Band %d...  occ = %lf vs. occ_min = %lf...\n",band+1,occ[band],occ_min);
                for (pw = 0; pw < npw; pw++) {
                    fread( & cg, sizeof(double), 1, f2);
                    //	          wfkin->cg[pw][band][0]=cg;
                    cg_band[pw][band][0] = cg;
                    fread( & cg, sizeof(double), 1, f2);
                    //	          wfkin->cg[pw][band][1]=cg;
                    cg_band[pw][band][1] = cg;
                }
                if (nspinor == 2) {
                    for (pw = 0; pw < npw; pw++) {
                        fread( & cg, sizeof(double), 1, f2);
                        cg_beta[pw][band][0] = cg;
                        fread( & cg, sizeof(double), 1, f2);
                        cg_beta[pw][band][1] = cg;
                    }
                }
                normalization = 0.0;
                if (occ[band] > occ_min) {
                    for (jx = 0; jx < den1.NGX - 1; jx++) {
                        for (jy = 0; jy < den1.NGY - 1; jy++) {
                            for (jz = 0; jz < den1.NGZ - 1; jz++) {
                                offsetRE = 2 * (jz + ngfftz * (jy + jx * ngffty));
                                offsetIM = (2 * (jz + ngfftz * (jy + jx * ngffty)) + 1);
                                *(psi_G_ptr + offsetRE) = (double) 0.0;
                                *(psi_G_ptr + offsetIM) = (double) 0.0;
                            }
                        }
                    }
                    for (pw1 = 0; pw1 < npw; pw1++) {
                        k_ngx = wfkin -> kg[pw1][0];
                        k_ngy = wfkin -> kg[pw1][1];
                        k_ngz = wfkin -> kg[pw1][2];
                        if (k_ngx < 0) k_ngx += ngfftx;
                        if (k_ngy < 0) k_ngy += ngffty;
                        if (k_ngz < 0) k_ngz += ngfftz;
                        offsetRE = 2 * (k_ngz + ngfftz * (k_ngy + k_ngx * ngffty));
                        offsetIM = (2 * (k_ngz + ngfftz * (k_ngy + k_ngx * ngffty)) + 1);
                        *(psi_G_ptr + offsetRE) = (double) cg_band[pw1][band][0];
                        *(psi_G_ptr + offsetIM) = (double) - cg_band[pw1][band][1];
                    } /* END PW1 LOOP */
                    fftw_execute(fftplan_psi);
                    //                 printf("   FFT done!\n");
                    for (atom = 0; atom < natom; atom++) {
                        for (l = 0; l < 4; l++) {
                            for (m = 0; m < 8; m++) {
                                for (j = 0; j < 3; j++) {
                                    Int_PjYlm_RE[l][m][j][atom] = 0.0;
                                    Int_PjYlm_IM[l][m][j][atom] = 0.0;
                                }
                            }
                        }
                    }
                    /* Make separable projections */
                    for (jz = 0; jz < ngz; jz++) {
                        /* fractional coordinates */
                        zf = (double) jz / (double) ngz;
                        for (jy = 0; jy < ngy; jy++) {
                            yf = (double) jy / (double) ngy;
                            for (jx = 0; jx < ngx; jx++) {
                                xf = (double) jx / (double) ngx;
                                voxcenter_x = xf * den1.cella_x + yf * den1.cellb_x + zf * den1.cellc_x;
                                voxcenter_y = xf * den1.cella_y + yf * den1.cellb_y + zf * den1.cellc_y;
                                voxcenter_z = xf * den1.cella_z + yf * den1.cellb_z + zf * den1.cellc_z;
                                offsetRE = 2 * (jz + ngfftz * (jy + jx * ngffty));
                                offsetIM = (2 * (jz + ngfftz * (jy + jx * ngffty)) + 1);
                                c1 = gsl_complex_rect((double) * (psi_R_ptr + offsetRE), (double) * (psi_R_ptr + offsetIM));
                                //                          c2 = gsl_complex_polar(1.0,-2*PI*(kpt[0][kptno]*jx/ngfftx+kpt[1][kptno]*jy/ngffty+kpt[2][kptno]*jz/ngfftz));
                                for (i = 0; i < vmap.neighcount2[jx][jy][jz]; i++) {
                                    ka = (vmap.ionmap[i][jx][jy][jz] >> 13 & 7) - 3;
                                    kb = (vmap.ionmap[i][jx][jy][jz] >> 10 & 7) - 3;
                                    kc = (vmap.ionmap[i][jx][jy][jz] >> 7 & 7) - 3;
                                    atom = vmap.ionmap[i][jx][jy][jz] & 127;
                                    xc = den1.Xcart[atom] + ka * den1.cella_x + kb * den1.cellb_x + kc * den1.cellc_x;
                                    yc = den1.Ycart[atom] + ka * den1.cella_y + kb * den1.cellb_y + kc * den1.cellc_y;
                                    zc = den1.Zcart[atom] + ka * den1.cella_z + kb * den1.cellb_z + kc * den1.cellc_z;
                                    c2 = gsl_complex_polar(1.0, -2 * PI * (kpt[0][kptno] * (1.0 * jx / ngfftx - 1.0 * ka) + kpt[1][kptno] * (1.0 * jy / ngffty - 1.0 * kb) + kpt[2][kptno] * (1.0 * jz / ngfftz - 1.0 * kc)));
                                    c1c2 = gsl_complex_mul(c1, c2);
                                    atom_type = typat[atom] - 1;
                                    xyz2sph(voxcenter_x - xc, voxcenter_y - yc, voxcenter_z - zc, & dist, & theta, & phi);
                                    costheta1 = cos(theta);
                                    for (m = 1; m < 4; m++) {
                                        sinmphi[m] = sin(m * phi);
                                        cosmphi[m] = cos(m * phi);
                                    }
                                    /* double Int_PjYlm[4][8][4][NATOMS_MAX];  /* l,m,j,atom */
                                    for (j = 0; j < 3; j++) {
                                        l = 0;
                                        Int_PjYlm_RE[0][0][j][atom] += vmap.Plj_RE[i][l][0][j][jx][jy][jz] * GSL_REAL(c1c2);
                                        Int_PjYlm_IM[0][0][j][atom] += vmap.Plj_RE[i][l][0][j][jx][jy][jz] * GSL_IMAG(c1c2);
                                        for (l = 1; l < psp.l_max[atom_type]; l++) {
                                            m = 0;
                                            Int_PjYlm_RE[l][m][j][atom] += vmap.Plj_RE[i][l][0][j][jx][jy][jz] * GSL_REAL(c1c2);
                                            Int_PjYlm_IM[l][m][j][atom] += vmap.Plj_RE[i][l][0][j][jx][jy][jz] * GSL_IMAG(c1c2);
                                            for (m = 1; m < l + 1; m++) {
                                                Int_PjYlm_RE[l][2 * m + 1][j][atom] += vmap.Plj_RE[i][l][2 * m + 1][j][jx][jy][jz] * GSL_REAL(c1c2);
                                                Int_PjYlm_RE[l][2 * m + 1][j][atom] += vmap.Plj_IM[i][l][2 * m + 1][j][jx][jy][jz] * GSL_IMAG(c1c2);
                                                Int_PjYlm_IM[l][2 * m + 1][j][atom] -= vmap.Plj_IM[i][l][2 * m + 1][j][jx][jy][jz] * GSL_REAL(c1c2);
                                                Int_PjYlm_IM[l][2 * m + 1][j][atom] += vmap.Plj_RE[i][l][2 * m + 1][j][jx][jy][jz] * GSL_IMAG(c1c2);
                                                Int_PjYlm_RE[l][2 * m][j][atom] += vmap.Plj_RE[i][l][2 * m][j][jx][jy][jz] * GSL_REAL(c1c2);
                                                Int_PjYlm_RE[l][2 * m][j][atom] += vmap.Plj_IM[i][l][2 * m][j][jx][jy][jz] * GSL_IMAG(c1c2);
                                                Int_PjYlm_IM[l][2 * m][j][atom] -= vmap.Plj_IM[i][l][2 * m][j][jx][jy][jz] * GSL_REAL(c1c2);
                                                Int_PjYlm_IM[l][2 * m][j][atom] += vmap.Plj_RE[i][l][2 * m][j][jx][jy][jz] * GSL_IMAG(c1c2);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    for (atom = 0; atom < natom; atom++) {
                        atom_type = typat[atom] - 1;
                        for (j = 0; j < 3; j++) {
                            for (i = 0; i < 3; i++) {
                                for (l = 0; l < psp.l_max[atom_type]; l++) {
                                    for (m = 0; m < 2 * l + 2; m++) {
                                        nonlocalE_byatom[dtset][atom][l] += occ[band] * Int_PjYlm_RE[l][m][j][atom] * pspin -> h[j][i][l][atom_type] * Int_PjYlm_RE[l][m][i][atom] * den1.VoxelV * den1.VoxelV / cellvolume;
                                        nonlocalE_byatom[dtset][atom][l] += occ[band] * Int_PjYlm_IM[l][m][j][atom] * pspin -> h[j][i][l][atom_type] * Int_PjYlm_IM[l][m][i][atom] * den1.VoxelV * den1.VoxelV / cellvolume;
                                    }
                                }
                            }
                        }
                    }
                } /* END IF OCC > 0 */
                fread( & j, sizeof(int), 1, f2);
            }
            for (atomno = 0; atomno < natom; atomno++) {
                nonlocalE_byatom[dtset][atomno][5] = 0.0;
                for (l = 0; l < 4; l++) {
                    nonlocalE_byatom[dtset][atomno][5] += nonlocalE_byatom[dtset][atomno][l];
                }
            }
            for (l = 0; l < natom; l++) {
                printf("    atom %d:  nonlocalE = %13.8lf = ", l + 1, nonlocalE_byatom[dtset][l][5]);
                printf(" (s) %13.8lf + ", nonlocalE_byatom[dtset][l][0]);
                printf(" (p) %13.8lf + ", nonlocalE_byatom[dtset][l][1]);
                printf(" (d) %13.8lf + ", nonlocalE_byatom[dtset][l][2]);
                printf(" (f) %13.8lf \n", nonlocalE_byatom[dtset][l][3]);
            }
            printf("\n");
        } /* END KPOINT LOOP */
    } /* END SPINPOL LOOP */
    fclose(f2);
    f3 = fopen(outfile, "a");
    for (l = 0; l < natom; l++) {
        printf("    atom %d:  nonlocalE = %13.8lf = ", l + 1, nonlocalE_byatom[dtset][l][5]);
        printf(" (s) %13.8lf + ", nonlocalE_byatom[dtset][l][0]);
        printf(" (p) %13.8lf + ", nonlocalE_byatom[dtset][l][1]);
        printf(" (d) %13.8lf + ", nonlocalE_byatom[dtset][l][2]);
        printf(" (f) %13.8lf \n", nonlocalE_byatom[dtset][l][3]);
        fprintf(f3, "  atom %d:  %20.17lf ", l + 1, nonlocalE_byatom[dtset][l][5]);
        fprintf(f3, " %20.17lf ", nonlocalE_byatom[dtset][l][0]);
        fprintf(f3, " %20.17lf ", nonlocalE_byatom[dtset][l][1]);
        fprintf(f3, " %20.17lf ", nonlocalE_byatom[dtset][l][2]);
        fprintf(f3, " %20.17lf \n", nonlocalE_byatom[dtset][l][3]);
    }
    fclose(f3);
}

void read_calc_den(char filename[200], struct XSFfile * XSFIN) {
    FILE * f2;
    int k;
    char codvsn[110];
    char title[132];
    char test;
    int headform;
    int fform;
    int bandtot;
    int date;
    int intxc;
    int ixc;
    int natom;
    int atomno, pw1, pw2;
    int ngfftx;
    int ngffty;
    int ngfftz;
    int nkpt;
    int npsp;
    int j = 0;
    int i, m, l;
    int jx, jy, jz;
    int stop = 0;
    int nspden;
    int nsppol;
    int ntypat;
    int occopt;
    int pertcase;
    int usepaw;
    double ecut;
    double ecutdg;
    double ecutsm;
    double ecut_eff;
    double qptnx;
    double qptny;
    double qptnz;
    double rprimd_ax;
    double rprimd_ay;
    double rprimd_az;
    double rprimd_bx;
    double rprimd_by;
    double rprimd_bz;
    double rprimd_cx;
    double rprimd_cy;
    double rprimd_cz;
    double ax_star;
    double ay_star;
    double az_star;
    double bx_star;
    double by_star;
    double bz_star;
    double cx_star;
    double cy_star;
    double cz_star;
    double ox, oy, oz;
    double theta1, phi1, g1;
    double theta2, phi2, g2;
    double stmbias;
    double tphysel;
    double tsmear;
    double znuclpsp;
    double zionpsp;
    double nonlocalE_byatom[NATOMS_MAX];
    double nonlocalE_byatom_IM[NATOMS_MAX];
    int pspso;
    int pspdat;
    int pspcod;
    int pspxc;
    int type;
    int lmn_size;
    int atom_type;
    int usewvl;
    int istwfk[KPOINTS_MAX];
    int istwfkv;
    int nband[BANDS_MAX];
    int nbandv;
    int npwarr[KPOINTS_MAX];
    int npwarrv;
    int so_psp[NATOMS_MAX];
    int symafm[MAX_SYM];
    double kpt[3][KPOINTS_MAX];
    double occ;
    double cgt[2][WAVES_MAX][2];
    double znucltypat[NATOMS_MAX];
    double wtk[KPOINTS_MAX];
    double wtkv;
    double residm;
    double x, y, z;
    double etotal, fermie;
    double Enonlocal_temp;
    double ga1, gb1, gc1;
    double gx1, gy1, gz1;
    double ga2, gb2, gc2;
    double gx2, gy2, gz2;
    double xred[3][NATOMS_MAX];
    double kptv;
    double php;
    int npw;
    int pw;
    int nspinor;
    int nband_temp;
    int kptno;
    int kx, ky, kz;
    int band;
    double normalization;
    double eigen, occ_temp, cg;
    double cellvolume, p1, p2;
    gsl_complex c1;
    gsl_complex c2;
    gsl_complex c1_star;
    gsl_complex c1c2;
    gsl_complex atom_phase;
    double atom_phase0;
    gsl_complex vg1g2;
    gsl_complex vg1g2_lm;
    gsl_complex vg1g2_temp;
    for (j = 0; j < NATOMS_MAX; j++) {
        nonlocalE_byatom[j] = 0.0;
        nonlocalE_byatom_IM[j] = 0.0;
    }
    f2 = fopen(filename, "rb+");
    if (f2 == NULL) {
        printf("%s not found.\n", filename);
        exit(0);
    }
    /*   READ HEADER OF WFK FILE    */
    j = 0;
    fread( & j, sizeof(int), 1, f2);
    j = fread(codvsn, sizeof(char), 6, f2);
    fread( & headform, sizeof(int), 1, f2);
    fread( & fform, sizeof(int), 1, f2);
    //    printf("%s %d %d\n",codvsn,headform,fform);
    fread( & j, sizeof(int), 1, f2);
    fread( & j, sizeof(int), 1, f2);
    fread( & bandtot, sizeof(int), 1, f2);
    /*    int date;*/
    fread( & date, sizeof(int), 1, f2);
    //    printf("%d %d\n",bandtot,date);
    /*    int intxc;  */
    fread( & intxc, sizeof(int), 1, f2);
    /*    int ixc;*/
    fread( & ixc, sizeof(int), 1, f2);
    /*    int natom;*/
    fread( & natom, sizeof(int), 1, f2);
    XSFIN -> NIONS = natom;
    /*    int ngfftx;*/
    fread( & ngfftx, sizeof(int), 1, f2);
    /*    int ngffty;*/
    fread( & ngffty, sizeof(int), 1, f2);
    /*    int ngfftz;*/
    fread( & ngfftz, sizeof(int), 1, f2);
    //    printf("ngfft = %d x %d x %d\n",ngfftx,ngffty,ngfftz);
    XSFIN -> NGX = ngfftx + 1;
    XSFIN -> NGY = ngffty + 1;
    XSFIN -> NGZ = ngfftz + 1;
    ngx = ngfftx;
    ngy = ngffty;
    ngz = ngfftz;
    AllocDbl(XSFIN);
    /*    int nkpt;*/
    fread( & nkpt, sizeof(int), 1, f2);
    /*    int nspden;*/
    fread( & nspden, sizeof(int), 1, f2);
    /*    int nspinor;*/
    fread( & nspinor, sizeof(int), 1, f2);
    /*    int nsppol;*/
    fread( & nsppol, sizeof(int), 1, f2);
    /*    int nsym;*/
    fread( & nsym, sizeof(int), 1, f2);
    /*    int npsp;*/
    fread( & npsp, sizeof(int), 1, f2);
    /*    int ntypat;*/
    fread( & ntypat, sizeof(int), 1, f2);
    /*    int occopt;*/
    fread( & occopt, sizeof(int), 1, f2);
    /*    int pertcase;*/
    fread( & pertcase, sizeof(int), 1, f2);
    /*    int usepaw;*/
    fread( & usepaw, sizeof(int), 1, f2);
    /*    double ecut;*/
    fread( & ecut, sizeof(double), 1, f2);
    /*     double ecutdg;*/
    fread( & ecutdg, sizeof(double), 1, f2);
    /*    double ecutsm;*/
    fread( & ecutsm, sizeof(double), 1, f2);
    /*    double ecut_eff;*/
    fread( & ecut_eff, sizeof(double), 1, f2);
    /*    double qptnx;*/
    fread( & qptnx, sizeof(double), 1, f2);
    /*    double qptny;*/
    fread( & qptny, sizeof(double), 1, f2);
    /*    double qptnz;*/
    fread( & qptnz, sizeof(double), 1, f2);
    /*    double rprimd_ax;*/
    fread( & rprimd_ax, sizeof(double), 1, f2);
    XSFIN -> cella_x = rprimd_ax;
    /*    double rprimd_ay;*/
    fread( & rprimd_ay, sizeof(double), 1, f2);
    XSFIN -> cella_y = rprimd_ay;
    /*    double rprimd_az;*/
    fread( & rprimd_az, sizeof(double), 1, f2);
    XSFIN -> cella_z = rprimd_az;
    /*    double rprimd_bx;*/
    fread( & rprimd_bx, sizeof(double), 1, f2);
    XSFIN -> cellb_x = rprimd_bx;
    /*    double rprimd_by;*/
    fread( & rprimd_by, sizeof(double), 1, f2);
    XSFIN -> cellb_y = rprimd_by;
    /*    double rprimd_bz;*/
    fread( & rprimd_bz, sizeof(double), 1, f2);
    XSFIN -> cellb_z = rprimd_bz;
    /*    double rprimd_cx;*/
    fread( & rprimd_cx, sizeof(double), 1, f2);
    XSFIN -> cellc_x = rprimd_cx;
    /*    double rprimd_cy;*/
    fread( & rprimd_cy, sizeof(double), 1, f2);
    XSFIN -> cellc_y = rprimd_cy;
    /*    double rprimd_cz;*/
    fread( & rprimd_cz, sizeof(double), 1, f2);
    XSFIN -> cellc_z = rprimd_cz;
    /*    double stmbias;*/
    cellvolume = (rprimd_ax * (rprimd_by * rprimd_cz - rprimd_bz * rprimd_cy) - rprimd_ay * (rprimd_bx * rprimd_cz - rprimd_bz * rprimd_cx) + rprimd_az * (rprimd_bx * rprimd_cy - rprimd_by * rprimd_cx));
    //    printf(" cella = %lf %lf %lf \n",rprimd_ax,rprimd_ay,rprimd_az);
    //    printf(" cellb = %lf %lf %lf \n",rprimd_bx,rprimd_by,rprimd_bz);
    //    printf(" cellc = %lf %lf %lf \n\n",rprimd_cx,rprimd_cy,rprimd_cz);
    //    printf(" cell volume = %lf\n\n",cellvolume);
    XSFIN -> cellvolume = cellvolume;
    XSFIN -> VoxelV = XSFIN -> cellvolume / ((XSFIN -> NGX - 1) * (XSFIN -> NGY - 1) * (XSFIN -> NGZ - 1));
    ax_star = 2 * PI * (rprimd_by * rprimd_cz - rprimd_cy * rprimd_bz) / cellvolume;
    ay_star = -2 * PI * (rprimd_bx * rprimd_cz - rprimd_cx * rprimd_bz) / cellvolume;
    az_star = 2 * PI * (rprimd_bx * rprimd_cy - rprimd_cx * rprimd_by) / cellvolume;
    bx_star = 2 * PI * (rprimd_cy * rprimd_az - rprimd_ay * rprimd_cz) / cellvolume;
    by_star = -2 * PI * (rprimd_cx * rprimd_az - rprimd_ax * rprimd_cz) / cellvolume;
    bz_star = 2 * PI * (rprimd_cx * rprimd_ay - rprimd_ax * rprimd_cy) / cellvolume;
    cx_star = 2 * PI * (rprimd_ay * rprimd_bz - rprimd_by * rprimd_az) / cellvolume;
    cy_star = -2 * PI * (rprimd_ax * rprimd_bz - rprimd_bx * rprimd_az) / cellvolume;
    cz_star = 2 * PI * (rprimd_ax * rprimd_by - rprimd_bx * rprimd_ay) / cellvolume;
    //    printf(" cella* = %lf %lf %lf \n",ax_star,ay_star,az_star);
    //    printf(" cellb* = %lf %lf %lf \n",bx_star,by_star,bz_star);
    //    printf(" cellc* = %lf %lf %lf \n\n",cx_star,cy_star,cz_star);
    //    printf(" a.a*=%lf, a.b*=%lf, a.c*=%lf\n",rprimd_ax*ax_star+rprimd_ay*ay_star+rprimd_az*az_star,rprimd_ax*bx_star+rprimd_ay*by_star+rprimd_az*bz_star,rprimd_ax*cx_star+rprimd_ay*cy_star+rprimd_az*cz_star);
    //    printf(" b.a*=%lf, b.b*=%lf, b.c*=%lf\n",rprimd_bx*ax_star+rprimd_by*ay_star+rprimd_bz*az_star,rprimd_bx*bx_star+rprimd_by*by_star+rprimd_bz*bz_star,rprimd_bx*cx_star+rprimd_by*cy_star+rprimd_bz*cz_star);
    //    printf(" c.a*=%lf, c.b*=%lf, c.c*=%lf\n\n",rprimd_cx*ax_star+rprimd_cy*ay_star+rprimd_cz*az_star,rprimd_cx*bx_star+rprimd_cy*by_star+rprimd_cz*bz_star,rprimd_cx*cx_star+rprimd_cy*cy_star+rprimd_cz*cz_star);
    fread( & stmbias, sizeof(double), 1, f2);
    /*    double tphysel;*/
    fread( & tphysel, sizeof(double), 1, f2);
    /*    double tsmear;*/
    fread( & tsmear, sizeof(double), 1, f2);
    /*    int usewvl;*/
    fread( & usewvl, sizeof(int), 1, f2);
    //    printf("natoms = %d   ecut = %lf  tsmear = %lf  occopt = %d \n",natom,ecut,tsmear, occopt);   
    fread( & j, sizeof(int), 1, f2);
    fread( & j, sizeof(int), 1, f2);
    /*    int istwfk [KPOINTS_MAX]; */
    //    printf("nkpt = %d    nsppol = %d   npsp = %d  ntypat = %d \n",nkpt,nsppol,npsp,ntypat);
    for (j = 0; j < nkpt; j++) {
        fread( & istwfkv, sizeof(int), 1, f2);
    }
    /*    int nband [BANDS_MAX]; */
    for (j = 0; j < (nkpt * nsppol); j++) {
        fread( & nbandv, sizeof(int), 1, f2);
    }
    /*    int npwarr [KPOINTS_MAX]; */
    for (j = 0; j < (nkpt); j++) {
        fread( & npwarrv, sizeof(int), 1, f2);
    }
    /*    int so_psp [NATOMS_MAX]; */
    for (j = 0; j < (npsp); j++) {
        fread( & so_psp[j], sizeof(int), 1, f2);
    }
    /*    int symafm [MAX_SYM]; */
    for (j = 0; j < (nsym); j++) {
        fread( & symafm[j], sizeof(int), 1, f2);
    }
    /*    int symrel [3][3][MAX_SYM]; */
    for (j = 0; j < (nsym); j++) {
        fread( & symrel[0][0][j], sizeof(int), 1, f2);
        fread( & symrel[1][0][j], sizeof(int), 1, f2);
        fread( & symrel[2][0][j], sizeof(int), 1, f2);
        fread( & symrel[0][1][j], sizeof(int), 1, f2);
        fread( & symrel[1][1][j], sizeof(int), 1, f2);
        fread( & symrel[2][1][j], sizeof(int), 1, f2);
        fread( & symrel[0][2][j], sizeof(int), 1, f2);
        fread( & symrel[1][2][j], sizeof(int), 1, f2);
        fread( & symrel[2][2][j], sizeof(int), 1, f2);
    }
    /*    int typat [NATOMS_MAX]; */
    for (j = 0; j < (natom); j++) {
        fread( & typat[j], sizeof(int), 1, f2);
        //	 printf("typeat = %d\n",typat[j]);         
    }
    /*    double kpt [3][KPOINTS_MAX]; */
    for (j = 0; j < (nkpt); j++) {
        fread( & kptv, sizeof(double), 1, f2);
        fread( & kptv, sizeof(double), 1, f2);
        fread( & kptv, sizeof(double), 1, f2);
        //         printf("kpoint %d:  %lf %lf %lf \n",j+1,kpt[0][j],kpt[1][j],kpt[2][j]);
    }
    /*    double occ(TBANDS_MAX); */
    for (j = 0; j < (bandtot); j++) {
        fread( & occ, sizeof(double), 1, f2);
    }
    /*    double tnons [3][MAX_SYM]; */
    for (j = 0; j < (nsym); j++) {
        fread( & tnons[0][j], sizeof(double), 1, f2);
        fread( & tnons[1][j], sizeof(double), 1, f2);
        fread( & tnons[2][j], sizeof(double), 1, f2);
    }
    /*    double znucltypat[NATOMS_MAX];  */
    for (j = 0; j < (ntypat); j++) {
        fread( & znucltypat[j], sizeof(double), 1, f2);
    }
    for (j = 0; j < natom; j++) {
        type = typat[j] - 1;
        //          printf("atom %d type=%d\n",j, type);
        XSFIN -> atomicno[j] = (int) znucltypat[type];
    }
    /*    double wtk[KPOINTS_MAX]; */
    for (j = 0; j < (nkpt); j++) {
        fread( & wtkv, sizeof(double), 1, f2);
    }
    fread( & j, sizeof(int), 1, f2);
    for (k = 0; k < (npsp); k++) {
        fread( & j, sizeof(int), 1, f2);
        fread(title, sizeof(char), 132, f2);
        //          printf("     %s\n",title);
        fread( & znuclpsp, sizeof(double), 1, f2);
        //          printf("     %lf ",znuclpsp);
        fread( & zionpsp, sizeof(double), 1, f2);
        //          printf("%lf ",zionpsp);
        fread( & pspso, sizeof(int), 1, f2);
        fread( & pspdat, sizeof(int), 1, f2);
        fread( & pspcod, sizeof(int), 1, f2);
        fread( & pspxc, sizeof(int), 1, f2);
        fread( & lmn_size, sizeof(int), 1, f2);
        //          printf("%d %d %d %d %d \n",pspso,pspdat,pspcod,pspxc,lmn_size);
        fread( & j, sizeof(int), 1, f2);
    }
    if (usepaw == 0) {
        fread( & j, sizeof(int), 1, f2);
        fread( & residm, sizeof(double), 1, f2);
        //        printf("     residm = %lf\n",residm);
        //          printf("Enter origin x y z: ");
        for (k = 0; k < natom; k++) {
            fread( & x, sizeof(double), 1, f2);
            fread( & y, sizeof(double), 1, f2);
            fread( & z, sizeof(double), 1, f2);
            //              printf("     Atom %d:  %lf %lf %lf \n",k+1,x,y,z);
            XSFIN -> Xcart[k] = (x) * XSFIN -> cella_x + (y) * XSFIN -> cellb_x + (z) * XSFIN -> cellc_x;
            XSFIN -> Ycart[k] = (x) * XSFIN -> cella_y + (y) * XSFIN -> cellb_y + (z) * XSFIN -> cellc_y;
            XSFIN -> Zcart[k] = (x) * XSFIN -> cella_z + (y) * XSFIN -> cellb_z + (z) * XSFIN -> cellc_z;
        }
        fread( & etotal, sizeof(double), 1, f2);
        fread( & fermie, sizeof(double), 1, f2);
        //          printf("     Etotal = %lf    FermiE = %lf\n\n",etotal,fermie);
        fread( & j, sizeof(int), 1, f2);
    } else {
        printf("Yikes!  usepaw!=0.  We haven't written code for this case yet. \n");
    }
    /*   HEADER FINISHED - READ DENSITY DATA   */
    fread( & j, sizeof(int), 1, f2);
    for (jz = 0; jz < XSFIN -> NGZ; jz++) {
        for (jy = 0; jy < XSFIN -> NGY; jy++) {
            for (jx = 0; jx < XSFIN -> NGX; jx++) {
                if ((jx < XSFIN -> NGX - 1) && (jy < XSFIN -> NGY - 1) && (jz < XSFIN -> NGZ - 1)) {
                    fread( & eigen, sizeof(double), 1, f2);
                    XSFIN -> grid[jx][jy][jz] = eigen;
                }
            }
        }
    }
    fread( & j, sizeof(int), 1, f2);
    fclose(f2);
    for (jz = 0; jz < XSFIN -> NGZ; jz++) {
        for (jy = 0; jy < XSFIN -> NGY; jy++) {
            for (jx = 0; jx < XSFIN -> NGX; jx++) {
                if (jx == XSFIN -> NGX - 1) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[0][jy][jz];
                if (jy == XSFIN -> NGY - 1) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[jx][0][jz];
                if (jz == XSFIN -> NGZ - 1) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[jx][jy][0];
                if ((jx == XSFIN -> NGX - 1) && (jy == XSFIN -> NGY - 1)) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[0][0][jz];
                if ((jx == XSFIN -> NGX - 1) && (jz == XSFIN -> NGZ - 1)) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[0][jy][0];
                if ((jy == XSFIN -> NGY - 1) && (jz == XSFIN -> NGZ - 1)) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[jx][0][0];

                if ((jx == XSFIN -> NGX - 1) && (jy == XSFIN -> NGY - 1) && (jz == XSFIN -> NGZ - 1)) XSFIN -> grid[jx][jy][jz] = XSFIN -> grid[0][0][0];

            }
        }
    }
}

void NL2XSF(struct XSFfile * XSFIN, struct XSFfile * XSFOUT1, struct XSFfile * XSFOUT2) {
    int i, jx, jy, jz, jx1, jy1, jz1, ha, hb, hc;
    int j, i1, i2, ix, iy, iz;
    double ix_d, iy_d, iz_d;
    int atom, atom_type;
    int l;
    int line_counter;
    int nvoxels_atm;
    int check;
    int stop;
    int stop2;
    int step;
    double r_max;
    double elementname[4];
    char outfilename[1000];
    char str[1000];
    int NGX, NGY, NGZ, k;
    double cell_volume, dist, dist2;
    double rcut, rcut2;
    double xf, yf, zf, xf2, yf2, zf2;
    double voxel_centerx;
    double voxel_centery;
    double voxel_centerz;
    double deltar;
    double Wnl_atom[NATOMS_MAX][4];
    double nelect_core;
    int n_epsatm_values = 0;
    int step1, step2;
    double averageE1_vs_r[10001];
    double averageE2_vs_r[10001];
    double sigmaE1_vs_r[10001];
    double sigmaE2_vs_r[10001];
    double NL1, NL2, NELECT;
    int averageP_vs_r_npoints[10001];
    int symm_error_flag = 0;
    int kptno;
    int do_map = 0;
    int atomic_no;
    double Z;
    double tolerance = 0.5 * pow(XSFIN -> VoxelV * pow(0.52917720859, 3.0), 0.33333333333333); //0.005-->0.2-->0.4-->0.5
    double Z_ion_temp;
    double E_core_up, E_core_down;
    int Zuse[100];
    double r_by_Z[100];
    FILE * f2;
    FILE * f4;
    char profile_file[200];
    double cella, cellb, cellc;
    NGX = XSFIN -> NGX;
    NGY = XSFIN -> NGY;
    NGZ = XSFIN -> NGZ;
    r_max = 4.0;
    cella = pow(XSFIN -> cella_x * XSFIN -> cella_x + XSFIN -> cella_y * XSFIN -> cella_y + XSFIN -> cella_z * XSFIN -> cella_z, 0.5) / (NGX - 1);
    cellb = pow(XSFIN -> cellb_x * XSFIN -> cellb_x + XSFIN -> cellb_y * XSFIN -> cellb_y + XSFIN -> cellb_z * XSFIN -> cellb_z, 0.5) / (NGY - 1);
    cellc = pow(XSFIN -> cellc_x * XSFIN -> cellc_x + XSFIN -> cellc_y * XSFIN -> cellc_y + XSFIN -> cellc_z * XSFIN -> cellc_z, 0.5) / (NGZ - 1);
    stop = 0;
    XSFOUT1 -> NGX = NGX;
    XSFOUT1 -> NGY = NGY;
    XSFOUT1 -> NGZ = NGZ;
    XSFOUT2 -> NGX = NGX;
    XSFOUT2 -> NGY = NGY;
    XSFOUT2 -> NGZ = NGZ;
    AllocDbl(XSFOUT1);
    AllocDbl(XSFOUT2);
    for (jz = 0; jz < NGZ; jz++) {
        for (jy = 0; jy < NGY; jy++) {
            for (jx = 0; jx < NGX; jx++) {
                XSFOUT1 -> grid[jx][jy][jz] = 0.0;
                XSFOUT2 -> grid[jx][jy][jz] = 0.0;
            }
        }
    }
    XSFOUT1 -> VoxelV = XSFIN -> VoxelV;
    XSFOUT2 -> VoxelV = XSFIN -> VoxelV;
    for (j = 0; j < XSFIN -> NIONS; j++) {
        Wnl_atom[j][0] = 0.0;
        Wnl_atom[j][1] = 0.0;
        Wnl_atom[j][2] = 0.0;
        Wnl_atom[j][3] = 0.0;
    }
    for (jz = 0; jz < NGZ - 1; jz++) {
        for (jy = 0; jy < NGY - 1; jy++) {
            for (jx = 0; jx < NGX - 1; jx++) {
                for (i = 0; i < vmap.neighcount2[jx][jy][jz]; i++) {
                    atom = vmap.ionmap[i][jx][jy][jz] & 127;
                    atom_type = typat[atom] - 1;
                    for (i1 = 0; i1 < 3; i1++) {
                        for (i2 = 0; i2 < 3; i2++) {
                            for (l = 0; l < l_max; l++) {
                                Wnl_atom[atom][l] += XSFIN -> grid[jx][jy][jz] * psp.h[i1][i2][l][atom_type] * vmap.Plj_RE[i][l][0][i1][jx][jy][jz] * vmap.Plj_RE[i][l][0][i2][jx][jy][jz] * XSFIN -> VoxelV;
                            }
                        }
                    }
                }
            } /*  end loop on jx */
        } /*  end loop on jy */
    } /*  end loop on jz */
    for (j = 0; j < XSFIN -> NIONS; j++) {
        printf("Weights:  (s) %lf (p) %lf (d) %lf (f) %lf\n", Wnl_atom[j][0], Wnl_atom[j][1], Wnl_atom[j][2], Wnl_atom[j][3]);
    }
    for (jz = 0; jz < NGZ - 1; jz++) {
        for (jy = 0; jy < NGY - 1; jy++) {
            for (jx = 0; jx < NGX - 1; jx++) {
                for (i = 0; i < vmap.neighcount2[jx][jy][jz]; i++) {
                    atom = vmap.ionmap[i][jx][jy][jz] & 127;
                    for (i1 = 0; i1 < 3; i1++) {
                        for (i2 = 0; i2 < 3; i2++) {
                            for (l = 0; l < l_max; l++) {
                                if (nonlocalE_byatom[0][atom][l] != 0.0) {
                                    XSFOUT2 -> grid[jx][jy][jz] += nonlocalE_byatom[0][atom][l] * XSFIN -> grid[jx][jy][jz] * psp.h[i1][i2][l][typat[atom] - 1] * vmap.Plj_RE[i][l][0][i1][jx][jy][jz] * vmap.Plj_RE[i][l][0][i2][jx][jy][jz] / Wnl_atom[atom][l];
                                }
                            }
                        }
                    }
                }
            } /*  end loop on jx */
        } /*  end loop on jy */
    } /*  end loop on jz */
    NELECT = 0.0;
    NL1 = 0.0;
    printf("\n");
    symm_error_flag = 0;
    printf("Applying space group symmetry operations...");
    for (jz = 0; jz < NGZ - 1; jz++) {
        zf = (double) jz / (double) ngz;
        for (jy = 0; jy < NGY - 1; jy++) {
            yf = (double) jy / (double) ngy;
            for (jx = 0; jx < NGX - 1; jx++) {
                xf = (double) jx / (double) ngx;
                for (j = 0; j < nsym; j++) {
                    xf2 = 1.0 * symrel[0][0][j] * xf + 1.0 * symrel[0][1][j] * yf + 1.0 * symrel[0][2][j] * zf + tnons[0][j];
                    yf2 = 1.0 * symrel[1][0][j] * xf + 1.0 * symrel[1][1][j] * yf + 1.0 * symrel[1][2][j] * zf + tnons[1][j];
                    zf2 = 1.0 * symrel[2][0][j] * xf + 1.0 * symrel[2][1][j] * yf + 1.0 * symrel[2][2][j] * zf + tnons[2][j];
                    ix_d = xf2 * (double) ngx;
                    iy_d = yf2 * (double) ngy;
                    iz_d = zf2 * (double) ngz;
                    stop = 0;
                    while (stop == 0) {
                        if (ix_d < -0.0000001) ix_d += 1.0 * ngx;
                        else stop = 1;
                    }
                    stop = 0;
                    while (stop == 0) {
                        if (iy_d < -0.0000001) iy_d += 1.0 * ngy;
                        else stop = 1;
                    }
                    stop = 0;
                    while (stop == 0) {
                        if (iz_d < -0.0000001) iz_d += 1.0 * ngz;
                        else stop = 1;
                    }
                    ix = (int)(floor(ix_d + 0.5));
                    iy = (int)(floor(iy_d + 0.5));
                    iz = (int)(floor(iz_d + 0.5));
                    if (fabs((double) ix - ix_d) > 0.001) symm_error_flag = 1;
                    if (fabs((double) iy - iy_d) > 0.001) symm_error_flag = 1;
                    if (fabs((double) iz - iz_d) > 0.001) symm_error_flag = 1;
                    ix = ix % ngx;
                    iy = iy % ngy;
                    iz = iz % ngz;
                    XSFOUT1 -> grid[ix][iy][iz] += XSFOUT2 -> grid[jx][jy][jz] / nsym;
                }
            }
        }
    }
    printf("Done!\n");
    for (jz = 0; jz < NGZ - 1; jz++) {
        for (jy = 0; jy < NGY - 1; jy++) {
            for (jx = 0; jx < NGX - 1; jx++) {
                NELECT += XSFIN -> grid[jx][jy][jz] * XSFIN -> VoxelV;
                NL1 += XSFOUT1 -> grid[jx][jy][jz] * XSFIN -> VoxelV;
            }
        }
    }
    printf("Total electron count in DS2_DEN file:  %lf\n", NELECT);
    printf("Total nonlocal energy for NL.xsf:  %lf\n", NL1);
    for (jy = 0; jy < XSFIN -> NGY - 1; jy++) {
        for (jx = 0; jx < XSFIN -> NGX - 1; jx++) {
            XSFOUT1 -> grid[jx][jy][NGZ - 1] = XSFOUT1 -> grid[jx][jy][0];
        }
    }
    for (jz = 0; jz < XSFIN -> NGZ - 1; jz++) {
        for (jx = 0; jx < XSFIN -> NGX - 1; jx++) {
            XSFOUT1 -> grid[jx][NGY - 1][jz] = XSFOUT1 -> grid[jx][0][jz];
        }
    }
    for (jz = 0; jz < XSFIN -> NGZ - 1; jz++) {
        for (jy = 0; jy < XSFIN -> NGY - 1; jy++) {
            XSFOUT1 -> grid[NGX - 1][jy][jz] = XSFOUT1 -> grid[0][jy][jz];
        }
    }
    for (jz = 0; jz < XSFIN -> NGZ - 1; jz++) {
        XSFOUT1 -> grid[NGX - 1][NGY - 1][jz] = XSFOUT1 -> grid[0][0][jz];
    }
    for (jy = 0; jy < XSFIN -> NGY - 1; jy++) {
        XSFOUT1 -> grid[NGX - 1][jy][NGZ - 1] = XSFOUT1 -> grid[0][jy][0];
    }
    for (jx = 0; jx < XSFIN -> NGX - 1; jx++) {
        XSFOUT1 -> grid[jx][NGY - 1][NGZ - 1] = XSFOUT1 -> grid[jx][0][0];
    }
    XSFOUT1 -> grid[NGX - 1][NGY - 1][NGZ - 1] = XSFOUT1 -> grid[0][0][0];
}

int main(int argc, char * argv[]) {
    int j;
    int nparams;
    FILE * f2;
    FILE * f4;
    double dV;
    int choice;
    char tim0[100];
    char tim1[100];
    char tim2[100];
    char KDENfile1[100];
    char KDENfile2[100];
    char POTfile1[100];
    char POTfile2[100];
    char DENfile0[100];
    char DENfile1[100];
    char DENfile2[100];
    char VHXCfile1[100];
    char VHXCfile2[100];
    char VHAfile1[100];
    char VHAfile2[100];
    char systemcmd[200];
    char WFKfilename[100];
    char DENfilename[100];
    char WFKfilename1[100];
    char WFKfilename2[100];
    char NLfilename1[100];
    char NLfilename2[100];
    char outfilename[100];
    char tmp1[100], tmp2[100];
    double res;
    double V1;
    double V2;
    int k_min, k_max;
    double Ealpha_map1, Ealpha_map2, Ealpha_nomap1, Ealpha_nomap2;
    double P_entropy;
    double P_Ewald;
    double P_nonlocal;
    double P_Ealpha;
    double min_occup;
    int dtsetnum;
    if (argc > 4) {
        strcpy(outfilename, argv[1]);
        strcpy(WFKfilename, argv[2]);
        sscanf(argv[4], "%lf", & min_occup);
        sscanf(argv[3], "%d", & dtsetnum);
    } else {
        printf("Usage:  nonlocal17 <_out file> <_o base> <dtsetnum> <minimum band occupancy> [optional: k-point]\n");
        exit(0);
    }
    if (argc > 5) {
        sscanf(argv[5], "%d", & k_use);
        printf("focusing on k-point %d\n", k_use);
    }
    read_psp_data(outfilename, & psp);
    sprintf(WFKfilename1, "%s_o_DS%d_WFK", WFKfilename, dtsetnum);
    sprintf(NLfilename1, "%s_o_DS%d_NL.xsf", WFKfilename, dtsetnum);
    strcpy(DENfilename, WFKfilename);
    if (dtsetnum == 1) strcat(DENfilename, "_o_DS1_DEN");
    if (dtsetnum == 2) strcat(DENfilename, "_o_DS2_DEN");
    if (dtsetnum == 3) strcat(DENfilename, "_o_DS3_DEN");
    read_calc_den(DENfilename, & den1);
    AllocInt( & vmap);
    Getkm( & den1);
    if (k_use > 0) sprintf(outfilename, "%s_o_DS%d_NL_kpt%d.log", outfilename, dtsetnum, k_use);
    else sprintf(outfilename, "%s_o_DS%d_NL.log", outfilename, dtsetnum);
    f2 = fopen(outfilename, "r");
    if (f2 == NULL) {
        f2 = fopen(outfilename, "w");
        fprintf(f2, "Nonlocal psp contribution to total energy from %s\n", WFKfilename1);
        fclose(f2);
        CoordSearchSphere( & den1, & vmap);
        read_calc_nonlocal(WFKfilename1, outfilename, & wfk1, & psp, min_occup, 0);
        strcpy(DENfilename, WFKfilename);
        strcat(DENfilename, "_o_DS2_DEN");
        read_calc_den(DENfilename, & den0);
    } else {
        printf("Found NL data file %s.\n", outfilename);
        fclose(f2);
        strcpy(DENfilename, WFKfilename);
        strcat(DENfilename, "_o_DS2_DEN");
        read_calc_den(DENfilename, & den0);
        printf("Reading nonlocal components...\n");
        f2 = fopen(outfilename, "r");
        finish_line(f2);
        for (j = 0; j < den0.NIONS; j++) {
            fscanf(f2, "%s %s %lf %lf %lf %lf %lf", tmp1, tmp2, & nonlocalE_byatom[0][j][5], & nonlocalE_byatom[0][j][0], & nonlocalE_byatom[0][j][1], & nonlocalE_byatom[0][j][2], & nonlocalE_byatom[0][j][3]);
            finish_line(f2);
            printf("    atom %d:  nonlocalE = %13.8lf = ", j + 1, nonlocalE_byatom[0][j][5]);
            printf(" (s) %13.8lf + ", nonlocalE_byatom[0][j][0]);
            printf(" (p) %13.8lf + ", nonlocalE_byatom[0][j][1]);
            printf(" (d) %13.8lf + ", nonlocalE_byatom[0][j][2]);
            printf(" (f) %13.8lf \n", nonlocalE_byatom[0][j][3]);
        }
        fclose(f2);
    }
    printf("Proceeding to map generation...\n");
    CoordSearchSphere2( & den0, & vmap);
    NL2XSF( & den0, & nlden1, & nlden2);
    f4 = fopen(NLfilename1, "w");
    outputXSF( & den0, & nlden1, f4);
    fclose(f4);
}
