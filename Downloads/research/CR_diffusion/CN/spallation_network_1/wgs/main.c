#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffusion_func.h"
int main(void) {

    double Rg[6] = {10., 16., 25., 40., 63., 100.}; // Rigidity
    double delta[6] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8}; // diffusion index
    double D0[10] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}; // log10 of diffusion coefficient
    double Rd[10] = {2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0}; // disk radius
    double ZH[6] = {3.0, 4.0, 5.0, 6.0, 7.0, 8.0}; // halo height
    double BC_ratio[6] = {0.2857, 0.2543, 0.219, 0.1834, 0.1589, 0.1247}; // B/C ratio at different rigidities
    double BC_ratio_err[6] = {0.0060531, 0.00529528, 0.00472969, 0.00430116, 0.00429535, 0.00410488};
    double BeC_ratio[6] = {0.0978, 0.088, 0.0769, 0.068, 0.0584, 0.0471}; // Be/C ratio at different rigidities
    double BeC_ratio_err[6] = {0.00233452, 0.00208806, 0.00189737, 0.00187883, 0.00180278, 0.0019105};
    initialize_grids_log();
    // Dynamically allocate memory for ndis and temp
    double *ndis_O16 = malloc(NR * NZ * sizeof(double));
    double *ndis_N15 = malloc(NR * NZ * sizeof(double));
    double *ndis_N14 = malloc(NR * NZ * sizeof(double));
    double *ndis_C13 = malloc(NR * NZ * sizeof(double));
    double *ndis_C12 = malloc(NR * NZ * sizeof(double));
    double *ndis_B11 = malloc(NR * NZ * sizeof(double));
    double *ndis_B10 = malloc(NR * NZ * sizeof(double));
    double *ndis_Be10 = malloc(NR * NZ * sizeof(double));
    double *ndis_Be9 = malloc(NR * NZ * sizeof(double));
    double *ndis_Be7 = malloc(NR * NZ * sizeof(double));
    double *ndis_H = malloc(NR * NZ * sizeof(double));
    double gas_type = -1.0, Dd;
    initialize_H(ndis_H, gas_type);

    char fn[100];
    snprintf(fn, sizeof(fn), "D%.2f_delta%.2f.bin",D0[0], delta[0]);

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 10; j++) {
            double chi2 = 0.0;
            for (int k = 0; k < 6; k++) {
                Dd = pow(10.0, 28.0 + D0[0]) * pow(Rg[k], delta[0]) / 3.086e21 / 3.086e21 * 3.15576e7; // in kpc^2 / yr;
                // Initialize particle density and apply boundary conditions, C and O are primary CRs, N and B are secondary CRs
                initialize_to_disk(ndis_O16, 1.0, dT, Rd[j]); // O16 are all primary
                initialize_to_disk(ndis_N14, 0.1, dT, Rd[j]); 
                initialize_to_disk(ndis_C12, 1.0, dT, Rd[j]); 
                initialize_to_0(ndis_N15);
                initialize_to_0(ndis_C13);
                initialize_to_0(ndis_B11);
                initialize_to_0(ndis_B10);
                initialize_to_0(ndis_Be10);
                initialize_to_0(ndis_Be9);
                initialize_to_0(ndis_Be7);
                // printf("Thread %d processing Rigidity: %f GV with Dd: %e\n", omp_get_thread_num(), Rg[i], Dd);        
                solve_diffusion_equation_CN(ndis_O16, ndis_N15, ndis_N14, ndis_C13, ndis_C12, ndis_B11, ndis_B10, ndis_Be10, ndis_Be9, ndis_Be7, ndis_H, Dd, dT, Rg[k], Rd[j], ZH[i]);
                double BC_ratio_M = (ndis_B11[80 * NZ] + ndis_B10[80 * NZ]) / (ndis_C13[80 * NZ] + ndis_C12[80 * NZ]);
                double BeC_ratio_M = (ndis_Be10[80 * NZ] + ndis_Be9[80 * NZ] + ndis_Be7[80 * NZ]) / (ndis_C13[80 * NZ] + ndis_C12[80 * NZ]);
                double BC_diff_M = BC_ratio_M - BC_ratio[k];
                double BeC_diff_M = BeC_ratio_M - BeC_ratio[k];
                chi2 += BC_diff_M * BC_diff_M / (BC_ratio_err[k] * BC_ratio_err[k]) + BeC_diff_M * BeC_diff_M / (BeC_ratio_err[k] * BeC_ratio_err[k]);
            }
            write_array_to_bin(fn, &chi2, 1);
        }
    }
    free(ndis_O16);
    free(ndis_N15);
    free(ndis_N14);
    free(ndis_C13);
    free(ndis_C12);
    free(ndis_B11);
    free(ndis_B10);
    free(ndis_Be10);
    free(ndis_Be9);
    free(ndis_Be7);
    free(ndis_H);
    return 0;
}



