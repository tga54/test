#include <stdlib.h>
#include <math.h>
#include "diffusion_func.h"
#include <string.h>  // for strncpy

const int NR1 = 100, NZ1 = 50;
const double dR1 = 0.1, dZ1 = 0.02, dR2 = 0.5, dZ2 = 0.1, Z_disk = 0.1;

typedef enum {
    O16 = 0, 
    N15 = 1,
    N14 = 2,
    C13 = 3,
    C12 = 4,
    B11 = 5,
    B10 = 6,
    Be10 = 7,
    Be9 = 8,
    Be7 = 9
} Isotope;

// cross section matrix; diagonal: total inelastic cross section; off-diagnal: spallation cross section from dragon2 
static double sigma[10][10] = {
//   O16         N15        N14         C13          C12        B11        B10       Be10       Be9        Be7
    {293.2e-27, 60.77e-27, 30.07e-27, 21.66e-27,  33.16e-27, 26.83e-27, 10.40e-27, 4.07e-27,  3.35e-27,  8.70e-27},  // O16
    {0.0,       279.7e-27, 20.15e-27, 46.14e-27,  43.82e-27, 36.48e-27, 13.59e-27, 3.40e-27,  7.98e-27,  5.34e-27},   // N15
    {0.0,       0.0,       265.9e-27, 16.26e-27,  52.00e-27, 29.22e-27, 10.16e-27, 2.70e-27,  2.71e-27,  12.1e-27},   // N14
    {0.0,       0.0,       0.0,       255.7e-27,  56.73e-27, 41.65e-27, 4.70e-27,  6.55e-27,  7.52e-27,  4.70e-27},   // C13
    {0.0,       0.0,       0.0,       0.0,        237.1e-27, 52.41e-27, 18.08e-27, 4.05e-27,  7.32e-27,  9.48e-27},   // C12
    {0.0,       0.0,       0.0,       0.0,        0.0,       211.9e-27, 40.85e-27, 12.95e-27, 15.15e-27, 4.56e-27},  // B11
    {0.0,       0.0,       0.0,       0.0,        0.0,       0.0,       197.1e-27, 0.0,       13.95e-27, 21.27e-27},  // B10
    {0.0,       0.0,       0.0,       0.0,        0.0,       0.0,       0.0,       197.0e-27, 36.37e-27, 36.94e-27},  // Be10
    {0.0,       0.0,       0.0,       0.0,        0.0,       0.0,       0.0,       0.0,       181.7e-27, 45.17e-27},   // Be9
    {0.0,       0.0,       0.0,       0.0,        0.0,       0.0,       0.0,       0.0,       0.0,       150.0e-27}   // Be7

};

static double Z_O = 8.0, Z_N = 7.0, Z_C = 6.0, Z_B = 5.0, Z_Be = 4.0, A_O16 = 16.0, A_N15 = 15.0, A_N14 = 14.0, A_C13 = 13.0, A_C12 = 12.0, A_B11 = 11.0, A_B10 = 10.0, A_Be10 = 10.0, A_Be9 = 9.0, A_Be7 = 7.0; // atomic number and mass number

static double c = 3e10, yr2s = 3.1536e7, ratio = 40.46285, tau_Be10 = 1.387e6, tau_stable = 1e100, rest_energy_p = 0.9315, enhancement_factor = 1.84; // tau in seconds, tau_stable for stable isotopes, rest_energy_p in GeV, enhancement_factor for cross-section to account for heavier ISM nuclei

double R[NR], Z[NZ], dR[NR], dZ[NZ];

void initialize_grids_log(void){
    double Rc = (double)NR1 * dR1;
    double Zc = (double)NZ1 * dZ1;
    double indR = log10(Rc);
    double indZ = log10(Zc);
    
    for (int i = 0; i < NR; i++) {
        if (i * dR1 <= Rc)
        {
            R[i] = i * dR1;
            dR[i] = dR1;
        } 
        else {
            R[i] = pow(10.0, ((double)i - (int)NR1) / 200.0 + indR ); // 400 steps for every order of magnitude, i.e. R[i] = 10**1.0000, R[i+1] = 10**1.0025, R[i+2] = 10**1.005, R3 = 10**1.0075....... 
            dR[i] = R[i] - R[i - 1]; 
        }
        // printf("R is %.5f and dR is %.5f : ", R[i], dR[i]);
    }
    for (int k = 0; k < NZ; k++) {
        if (k * dZ1 <= Zc){
            Z[k] = k * dZ1;
            dZ[k] = dZ1;
        } 
        else{
            Z[k] = pow(10.0, ((double)k - (int)NZ1) / 50.0 + indZ ); //  100 steps for every order of magnitude: i.e. R[i] = 10**1.00, R[i+1] = 10**1.01, R[i+2] = 10**1.02, R[i+3] = 10**1.03....... 
            dZ[k] = Z[k] - Z[k - 1];
        }
        // printf("Z is %.5f and dZ is %.5f : ", Z[k], dZ[k]);
    }
}

void initialize_grids_linear(void){
    double Rc = (double)NR1 * dR1;
    double Zc = (double)NZ1 * dZ1;
    for (int i = 0; i < NR; i++) {
        if (i * dR1 <= Rc)
        {
            R[i] = i * dR1;
            dR[i] = dR1;
        } 
        else {
            R[i] = Rc + (i - NR1) * dR2; 
            dR[i] = dR2; 
        }
    }
    for (int k = 0; k < NZ; k++) {
        if (k * dZ1 <= Zc){
            Z[k] = k * dZ1;
            dZ[k] = dZ1;
        } 
        else{
            Z[k] = Zc + (k - NZ1) * dZ2;  
            dZ[k] = dZ2;
        }
    }
}

void initialize_to_disk(double *ndis, const double norm, const double dT, const double R_disk) {
    for (int i = 0; i < NR; i++)
        for (int j = 0; j < NZ; j++)
            ndis[i * NZ + j] = (R[i] <= R_disk && Z[j] <= Z_disk) ? norm * dT : 0.0;
}

void initialize_to_src(double *ndis, const double dT) {
    for (int i = 0; i < NR; i++){
        for (int j = 0; j < NZ; j++){
	    if (Z[j] < Z_disk){
	        ndis[i * NZ + j] = 64.6 * pow(R[i], 2.35) * exp(- R[i] / 1.528 ) * dT;
	    }
        } 
    }
}

void initialize_to_0(double *ndis) {
    for (int i = 0; i < NR; i++)
        for (int j = 0; j < NZ; j++)
            ndis[i * NZ + j] = 0.0;
}

void initialize_from_file(const char *filename, double *ndis) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening input file");
        return;
    }

    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            if (fscanf(file, "%lf,", &ndis[i * NZ + j]) != 1) {
                fprintf(stderr, "Error reading element [%d][%d]\n", i, j);
                fclose(file);
                return;
            }
        }
        // Consume the newline character after each row
        int c = fgetc(file);
        if (c != '\n' && c != EOF) {
            ungetc(c, file);  // put back if not newline or EOF
        }
    }
    fclose(file);
}


void apply_boundary_conditions(double *ndis, double ZH) {

    // --- R = 0 (reflective) ---
    for (int k = 0; k < NZ; k++) {
        ndis[0 * NZ + k] = ndis[1 * NZ + k];
    }

    // --- R = Rmax (absorbing) ---
    for (int k = 0; k < NZ; k++) {
        ndis[(NR - 1) * NZ + k] = 0.0;
    }

    // --- Z = 0 (reflective) ---
    for (int j = 0; j < NR; j++) {
        ndis[j * NZ + 0] = ndis[j * NZ + 1];
    }

    // --- Z = ZH (absorbing boundary) ---
    for (int j = 0; j < NR; j++) {
        for (int k = 0; k < NZ; k++) {
            if (Z[k] >= ZH) {
                ndis[j * NZ + k] = 0.0;
            }
        }
    }
}

void initialize_H(double *ndis_H, const double nH) {  
// nH = 0.0: no spallation; nH > 0: flat hydrogen distribution; nH = -1.0:we take the distribution of H_2 and HI from doi:10.1093/mnras/staa1017, and HII from doi: 10.1111/j.1365-2966.2004.08349.x 
    if (nH >= 0.0){
        for (int i = 0; i < NR; i++) {
            for (int j = 0; j < NZ; j++) {
                ndis_H[i*NZ+j] = nH;
            }
        }
    }
    else if (nH == -1.0){
        double *ndis_H2 = malloc(NR * NZ * sizeof(double));
        double *ndis_HI = malloc(NR * NZ * sizeof(double));
        double *ndis_HII = malloc(NR * NZ * sizeof(double));
        for (int i = 0; i < NR; i++) {
            for (int j = 0; j < NZ; j++) {
                double r0 = sqrt(R[i]*R[i] + Z[j]*Z[j]);
                ndis_H2[i*NZ+j] = ratio * 2200. / 4. / 45. * exp(-12./R[i] - R[i] / 1.5) / pow(cosh(Z[j] / 2. / 0.045), 2);
                ndis_HI[i*NZ+j] = ratio * 53. / 4. / 85. * exp(-4./R[i] - R[i] / 7) / pow(cosh(Z[j] / 2. / 0.085), 2);
                ndis_HII[i*NZ+j] = 0.00015 * (1. + 3.7*log(1 + r0 / 20.) / (r0/20.) - 1.0277);
                ndis_H[i*NZ+j] =  2.0 * ndis_H2[i*NZ+j] + ndis_HI[i*NZ+j] + ndis_HII[i*NZ+j];
                // ndis_H[i*NZ+j] = ndis_HII[i*NZ+j];
            }
        }
        free(ndis_H2); free(ndis_HI); free(ndis_HII);
    }
    else{
        perror("Error input of nH!");
        return;
    }
}


void write_array_to_bin(const char *filename, const double *arr, const int size) {
    FILE *file = fopen(filename, "ab");
    if (!file) {
        perror("File opening failed");
        return;
    }
    fwrite(arr, sizeof(double), size, file);
    fclose(file);
}

static inline void thomas_solve(const double *a, const double *b, const double *c, const double *d, double *x, double *cp, double *dp, int n)
{
    if (n < 2) {
        if (n == 1) x[0] = d[0] / b[0];
        return;
    }
    // Modify the first coefficients
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    // Forward sweep
    for (int i = 1; i < n; ++i) {
        double inv_denom = 1.0 / (b[i] - a[i] * cp[i - 1]);
        // if (denom == 0.0) {
        //     fprintf(stderr, "Zero pivot in Thomas algorithm at i=%d.\n", i);
        //     free(cp); free(dp);
        //     return;
        // }
        cp[i] = (i < n - 1) ? c[i] * inv_denom : 0.0;
        dp[i] = (d[i] - a[i] * dp[i - 1]) * inv_denom;
    }
    // Back substitution
    x[n - 1] = dp[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
}

void CN_scheme(double *ndis, double *src, double *ndis_H, double *dT_tau, double *temp, double beta, double sigma, double tau_decay, double D, double dT, double ZH) {

    // define ADI coefficients
    double alpha_R[NR], alpha_Z[NZ];
    for (int i = 0; i < NR - 1; i++) alpha_R[i] = D * dT / (dR[i] + dR[i + 1]) ;
    for (int i = 0; i < NZ - 1; i++) alpha_Z[i] = D * dT / (dZ[i] + dZ[i + 1]) ;

    // fragmentation and decay term
    for (int idx = 0; idx < NR * NZ; idx++) {            
        dT_tau[idx] = dT * yr2s * (ndis_H[idx] * sigma * beta * c + 1.0 / tau_decay); 
    } 

    // Step 1: Implicit R, Explicit Z
    // Thread-private arrays (allocated on stack, ~3.5KB per thread, safe)
    double CN_Ra[NR], CN_Rb[NR], CN_Rc[NR], CN_Rd[NR], CN_Rx[NR], CN_Rcp[NR], CN_Rdp[NR];
    for (int k = 1; k < NZ - 1; k++) {     
        for (int j = 1; j < NR - 1; j++) {
            CN_Rd[j] = alpha_Z[k] / dZ[k] * ndis[j * NZ + k - 1] + (1.0 - alpha_Z[k] * (dZ[k + 1] + dZ[k]) / dZ[k + 1] / dZ[k]) * ndis[j * NZ + k] + alpha_Z[k] / dZ[k + 1] * ndis[j * NZ + k + 1] - ndis[j * NZ + k] * dT_tau[j * NZ + k] / 4.0;
        }                   // right hand side
        for (int j = 1; j < NR - 1; j++){     // cylinderical coordinate system  
            CN_Ra[j] = - alpha_R[j] * (1.0 / dR[j] - 0.5 / R[j]);
            CN_Rb[j] = 1.0 + alpha_R[j] * (dR[j + 1] + dR[j]) / dR[j + 1] / dR[j] + dT_tau[j * NZ + k] / 4.0;
            CN_Rc[j] = - alpha_R[j] * (1.0 / dR[j + 1] + 0.5 / R[j]);
        }
        CN_Ra[0] = 0.0;
        CN_Rb[0] = 1.0;
        CN_Rc[0] = -1.0;
        CN_Rd[0] = 0.0;
        thomas_solve(CN_Ra, CN_Rb, CN_Rc, CN_Rd, CN_Rx, CN_Rcp, CN_Rdp, NR - 1);
        CN_Rx[NR - 1] = 0.0; // enforce absorbing boundary at Rmax
        for (int j = 0; j < NR; j++) temp[j * NZ + k] = CN_Rx[j];
    }
    apply_boundary_conditions(temp, ZH);
    // Step 2: Implicit Z, Explicit R
    double CN_Za[NZ], CN_Zb[NZ], CN_Zc[NZ], CN_Zd[NZ], CN_Zx[NZ], CN_Zcp[NZ], CN_Zdp[NZ];
    for (int j = 1; j < NR - 1; j++){
        for (int k = 1; k < NZ - 1; k++){
            CN_Zd[k] = alpha_R[j] * (1.0 / dR[j] - 0.5 / R[j]) * temp[(j-1) * NZ + k] + (1.0 - alpha_R[j] * (dR[j + 1] + dR[j]) / dR[j + 1] / dR[j]) * temp[j * NZ + k] + alpha_R[j] * (1.0 / dR[j + 1] + 0.5 / R[j]) * temp[(j+1) * NZ + k] - temp[j * NZ + k] * dT_tau[j * NZ + k] / 4.0;
        }
        for (int k = 1; k < NZ - 1; k++){
            CN_Za[k] = - alpha_Z[k] / dZ[k];
            CN_Zb[k] = 1.0 + alpha_Z[k] * (dZ[k + 1] + dZ[k]) / dZ[k + 1] / dZ[k] + dT_tau[j * NZ + k] / 4.0;
            CN_Zc[k] = - alpha_Z[k] / dZ[k + 1];
        }
        CN_Za[0] = 0.0;
        CN_Zb[0] = 1.0;
        CN_Zc[0] = -1.0;
        CN_Zd[0] = 0.0;
        thomas_solve(CN_Za, CN_Zb, CN_Zc, CN_Zd, CN_Zx, CN_Zcp, CN_Zdp, NZ - 1);
        CN_Zx[NZ - 1] = 0.0; // enforce absorbing boundary at ZH
        for (int k = 0; k < NZ; k++) ndis[j * NZ + k] = CN_Zx[k] + src[j * NZ + k];
    }
    apply_boundary_conditions(ndis, ZH);
}


void solve_diffusion_equation_CN(double *ndis_O16, double *ndis_N15, double *ndis_N14, double *ndis_C13, double *ndis_C12, double *ndis_B11, double *ndis_B10, double *ndis_Be10, double *ndis_Be9, double *ndis_Be7, double *ndis_H, double D, const double dT, double Rg, double R_disk, double ZH) {
    double *src_O16 = malloc(NR * NZ * sizeof(double));
    double *src_N15 = malloc(NR * NZ * sizeof(double));
    double *src_N14 = malloc(NR * NZ * sizeof(double));
    double *src_C13 = malloc(NR * NZ * sizeof(double));
    double *src_C12 = malloc(NR * NZ * sizeof(double));
    double *src_B11 = malloc(NR * NZ * sizeof(double));
    double *src_B10 = malloc(NR * NZ * sizeof(double));
    double *src_Be10 = malloc(NR * NZ * sizeof(double));
    double *src_Be9 = malloc(NR * NZ * sizeof(double));
    double *src_Be7 = malloc(NR * NZ * sizeof(double));
    double *temp = malloc(NR * NZ * sizeof(double));
    double *dT_tau = malloc(NR * NZ * sizeof(double));

    // O16 are all primary 
    initialize_to_disk(src_O16, 1.0, dT, R_disk); 

    // these are all secondary
    initialize_to_0(src_N15);
    initialize_to_0(src_C13);
    initialize_to_0(src_B11);
    initialize_to_0(src_B10);
    initialize_to_0(src_Be10);
    initialize_to_0(src_Be9);
    initialize_to_0(src_Be7);
    initialize_to_0(temp);
    initialize_to_0(dT_tau);

    double beta_O16 = Z_O * Rg / sqrt(Z_O * Z_O * Rg * Rg + A_O16 * A_O16 * rest_energy_p * rest_energy_p); // beta = v/c 
    double beta_N15 = Z_N * Rg / sqrt(Z_N * Z_N * Rg * Rg + A_N15 * A_N15 * rest_energy_p * rest_energy_p); // beta = v/c 
    double beta_N14 = Z_N * Rg / sqrt(Z_N * Z_N * Rg * Rg + A_N14 * A_N14 * rest_energy_p * rest_energy_p); // beta = v/c 
    double beta_C13 = Z_C * Rg / sqrt(Z_C * Z_C * Rg * Rg + A_C13 * A_C13 * rest_energy_p * rest_energy_p);
    double beta_C12 = Z_C * Rg / sqrt(Z_C * Z_C * Rg * Rg + A_C12 * A_C12 * rest_energy_p * rest_energy_p);
    double beta_B11 = Z_B * Rg / sqrt(Z_B * Z_B * Rg * Rg + A_B11 * A_B11 * rest_energy_p * rest_energy_p);
    double beta_B10 = Z_B * Rg / sqrt(Z_B * Z_B * Rg * Rg + A_B10 * A_B10 * rest_energy_p * rest_energy_p);
    double beta_Be10 = Z_Be * Rg / sqrt(Z_Be * Z_Be * Rg * Rg + A_Be10 * A_Be10 * rest_energy_p * rest_energy_p);
    double beta_Be9 = Z_Be * Rg / sqrt(Z_Be * Z_Be * Rg * Rg + A_Be9 * A_Be9 * rest_energy_p * rest_energy_p);
    double beta_Be7 = Z_Be * Rg / sqrt(Z_Be * Z_Be * Rg * Rg + A_Be7 * A_Be7 * rest_energy_p * rest_energy_p);

    for (int t = 0; t < NT; t++){
        // sigma is cross-section matrix defined above, sigma[a][a] is the total inelastic cross-section of a
        CN_scheme(ndis_O16, src_O16, ndis_H, dT_tau, temp, beta_O16, sigma[O16][O16]* enhancement_factor, tau_stable, beta_O16 * D, dT, ZH);
        // N14, C12 are combined of primary and secondary, suppose primary source terms are 10% and 100% of O16 respectively, refer to doi:10.1016/j.physrep.2020.09.003
        initialize_to_disk(src_N14, 0.1, dT, R_disk);
        initialize_to_disk(src_C12, 1.0, dT, R_disk);
        double l_year = c * dT * yr2s; // path length in one time step
        double gamma = sqrt(1.0 + pow(Rg * 4.0 / 10.0 / rest_energy_p, 2.0)); // Lorentz factor for Be-10 at rigidity Rg (GV)
        double tau_Be10_eff = tau_Be10 * yr2s * gamma; // time dilation, unit second

        for (int idx = 0; idx < NR * NZ; idx++){
            // spallation source terms, sigma is cross-section matrix defined above, sigma[a][b]: from a to b
            src_N15[idx] = ndis_O16[idx] * sigma[O16][N15] * ndis_H[idx] * beta_O16 * l_year * enhancement_factor;  
            src_N14[idx] += (ndis_O16[idx] * sigma[O16][N14] * beta_O16 + ndis_N15[idx] * sigma[N15][N14]* beta_N15) * ndis_H[idx] * l_year* enhancement_factor; 
            src_C13[idx] = (ndis_O16[idx] * sigma[O16][C13] * beta_O16 + ndis_N15[idx] * sigma[N15][C13] * beta_N15
                          + ndis_N14[idx] * sigma[N14][C13]* beta_N14) * ndis_H[idx] * l_year* enhancement_factor; 
            src_C12[idx] += (ndis_O16[idx] * sigma[O16][C12]* beta_O16 + ndis_N15[idx] * sigma[N15][C12] * beta_N15
                           + ndis_N14[idx] * sigma[N14][C12] * beta_N14 + ndis_C13[idx] * sigma[C13][C12] * beta_C13) * ndis_H[idx] * l_year* enhancement_factor;
            src_B11[idx] = (ndis_O16[idx] * sigma[O16][B11]* beta_O16 + ndis_N15[idx] * sigma[N15][B11] * beta_N15
                          + ndis_N14[idx] * sigma[N14][B11]* beta_N14 + ndis_C13[idx] * sigma[C13][B11] * beta_C13
                          + ndis_C12[idx] * sigma[C12][B11]* beta_C12) * ndis_H[idx] * l_year* enhancement_factor; 
            src_B10[idx] = (ndis_O16[idx] * sigma[O16][B10]* beta_O16 + ndis_N15[idx] * sigma[N15][B10] * beta_N15
                          + ndis_N14[idx] * sigma[N14][B10]* beta_N14 + ndis_C13[idx] * sigma[C13][B10] * beta_C13
                          + ndis_C12[idx] * sigma[C12][B10]* beta_C12 + ndis_B11[idx] * sigma[B11][B10] * beta_B11) * ndis_H[idx] * l_year* enhancement_factor + ndis_Be10[idx] / tau_Be10_eff * dT * yr2s;   
            src_Be10[idx] = (ndis_O16[idx] * sigma[O16][Be10]* beta_O16 + ndis_N15[idx] * sigma[N15][Be10] * beta_N15 
                           + ndis_N14[idx] * sigma[N14][Be10]* beta_N14 + ndis_C13[idx] * sigma[C13][Be10] * beta_C13 
                           + ndis_C12[idx] * sigma[C12][Be10]* beta_C12 + ndis_B11[idx] * sigma[B11][Be10] * beta_B11) * ndis_H[idx] * l_year* enhancement_factor; 

            src_Be9[idx] = (ndis_O16[idx] * sigma[O16][Be9] * beta_O16 + ndis_N15[idx] * sigma[N15][Be9]* beta_N15
                          + ndis_N14[idx] * sigma[N14][Be9] * beta_N14 + ndis_C13[idx] * sigma[C13][Be9]* beta_C13
                          + ndis_C12[idx] * sigma[C12][Be9] * beta_C12 + ndis_B11[idx] * sigma[B11][Be9]* beta_B11
                          + ndis_B10[idx] * sigma[B10][Be9] * beta_B10+ ndis_Be10[idx] * sigma[Be10][Be9] * beta_Be10) * ndis_H[idx] * l_year* enhancement_factor;
            src_Be7[idx] = (ndis_O16[idx] * sigma[O16][Be7] * beta_O16 + ndis_N15[idx] * sigma[N15][Be7]* beta_N15
                          + ndis_N14[idx] * sigma[N14][Be7] * beta_N14 + ndis_C13[idx] * sigma[C13][Be7]* beta_C13
                          + ndis_C12[idx] * sigma[C12][Be7] * beta_C12 + ndis_B11[idx] * sigma[B11][Be7]* beta_B11
                          + ndis_B10[idx] * sigma[B10][Be7] * beta_B10+ ndis_Be10[idx] * sigma[Be10][Be7] * beta_Be10
                          + ndis_Be9[idx] * sigma[Be9][Be7] * beta_Be9) * ndis_H[idx] * l_year* enhancement_factor;
        }


        // sigma is cross-section matrix defined above, sigma[a][a]: is the total inelastic cross-section of a
        CN_scheme(ndis_N15, src_N15, ndis_H, dT_tau, temp, beta_N15, sigma[N15][N15]* enhancement_factor, tau_stable, beta_N15 * D, dT, ZH);
        CN_scheme(ndis_N14, src_N14, ndis_H, dT_tau, temp, beta_N14, sigma[N14][N14]* enhancement_factor, tau_stable, beta_N14 * D, dT, ZH);
        CN_scheme(ndis_C13, src_C13, ndis_H, dT_tau, temp, beta_C13, sigma[C13][C13]* enhancement_factor, tau_stable, beta_C13 * D, dT, ZH);
        CN_scheme(ndis_C12, src_C12, ndis_H, dT_tau, temp, beta_C12, sigma[C12][C12]* enhancement_factor, tau_stable, beta_C12 * D, dT, ZH);
        CN_scheme(ndis_B11, src_B11, ndis_H, dT_tau, temp, beta_B11, sigma[B11][B11]* enhancement_factor, tau_stable, beta_B11 * D, dT, ZH);
        CN_scheme(ndis_B10, src_B10, ndis_H, dT_tau, temp, beta_B10, sigma[B10][B10]* enhancement_factor, tau_stable, beta_B10 * D, dT, ZH);
        CN_scheme(ndis_Be10, src_Be10, ndis_H, dT_tau, temp, beta_Be10, sigma[Be10][Be10]* enhancement_factor, tau_Be10_eff, beta_Be10 * D, dT, ZH);
        CN_scheme(ndis_Be9, src_Be9, ndis_H, dT_tau, temp, beta_Be9, sigma[Be9][Be9]* enhancement_factor, tau_stable, beta_Be9 * D, dT, ZH);
        CN_scheme(ndis_Be7, src_Be7, ndis_H, dT_tau, temp, beta_Be7, sigma[Be7][Be7]* enhancement_factor, tau_stable, beta_Be7 * D, dT, ZH);
    }
    free(temp);
    free(dT_tau);
    free(src_O16);
    free(src_N15);
    free(src_N14);
    free(src_C13);
    free(src_C12);
    free(src_B11);
    free(src_B10);
    free(src_Be10);
    free(src_Be9);
    free(src_Be7);
}

void write_array2D_to_txt(const char *filename, const double *arr) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening output file");
        return;
    }

    for (int i = 0; i < NR; i++) {
        for (int j = 0; j < NZ; j++) {
            fprintf(file, "%f,", arr[i * NZ + j]);
        }
        fputc('\n', file);  // new row
    }
    fclose(file);
}

void write_array1D_to_txt(const char *filename, const double *arr) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening output file");
        return;
    }
    // double ind_end = log10((double)NT) * 10.0;  

    for (int i = 0; i < NT; i++) {
        // double ind = pow(10.0, i / 10.0 + 1.0); //write 10 values for every order of magnitude, and start from N[10];
        fprintf(file, "%f,",  arr[i]);
    }
    fclose(file);
}

