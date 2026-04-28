#ifndef FUNC_H
#define FUNC_H

#include <stdio.h>

#define NR 161
#define NZ 101
#define NT 2001
#define dT 50000.0 
extern const int NR1, NZ1;
extern const double dR1, dZ1, dR2, dZ2, R_disk, Z_disk, R_inj;
extern double R[NR], Z[NZ], dR[NR], dZ[NZ];

typedef struct {
    int NR_BH, NZ_BH;
    double dR_BH, dZ_BH, dT_BH, Dd_BH, Dh_BH, ss_time_BH, ndis_H_BH, R_disk_BH;
    char description_BH[128];
} BinHeader;


void initialize_grids_linear(void);
void initialize_grids_log(void);
void initialize_to_disk(double *ndis, const double norm, const double dT, const double R_disk);
void initialize_to_src(double *ndis, const double dT);
void initialize_to_0(double *ndis);
void initialize_from_file(const char *filename, double *ndis);
void initialize_H(double *ndis_H, const double nH);
void apply_boundary_conditions(double *ndis, double ZH);

void CN_scheme(double *ndis, double *src, double *ndis_H, double *dT_tau, double *temp, double beta, double sigma, double tau_decay, double D, double dT);
void solve_diffusion_equation_CN(double *ndis_O16, double *ndis_N15, double *ndis_N14, double *ndis_C13, double *ndis_C12, double *ndis_B11, double *ndis_B10, double *ndis_Be10, double *ndis_Be9, double *ndis_Be7, double *ndis_H, double D, const double dT, double Rg, double R_disk, double ZH);
void write_array2D_to_txt(const char *filename, const double *ndis);
void write_array1D_to_txt(const char *filename, const double *N);
void write_array_to_bin(const char *filename, const double *arr, const int size);
void write_array_header_to_bin(const char *filename, const double *arr, const int size, const double dT, const double *ss_time, const double Dd, const double Dh, const double gas_type, const char *grid_type);

#endif
