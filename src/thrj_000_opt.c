/*
 * thrj_000_opt.c
 * 
 * Description: Optimized C code for calculating the coupling matrix for CMB temperature maps
 * Author: Georgia Kiddier
 * Date: 18-02-2025
 * License: MIT (or specify your chosen license)
 * 
 * Compilation:
 *     /opt/homebrew/bin/gcc-14 -O3 -march=native -o thrj_000_opt thrj_000_opt.c -lm
 * 
 * Usage:
 *     ./thrj_000_opt <power_spectrum_file> <output_matrix_file> <lmax> <write_output>
 * 
 * This program computes the coupling matrix M_l1l2 using precomputed 3j symbols.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

// Global variable to store the maximum multipole
int lmax;

/*
 * Function to calculate the 3j symbol (0 0 0) using precomputed arrays.
 */
double calculate_threej_000(int J, int J1minus, int J2minus, int J3minus, double* one_jp1, double* g, double* one_g) {
    //if (J % 2 == 0) { // 3j symbol is nonzero only if J is even

        int half_J1minus = J1minus / 2;
        int half_J2minus = J2minus / 2;
        int half_J3minus = J3minus / 2;
        int half_J = J / 2;

        double sign = (half_J % 2) ? -1.0 : 1.0;
        return sign * sqrt(one_jp1[J] * g[half_J1minus] * g[half_J2minus] * g[half_J3minus] * one_g[half_J]);
}


int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <input_file> <output_file> <lmax> <write_output (yes/no)>\n", argv[0]);
        return 1;
    }     
    lmax=atoi(argv[3]);
    char *write_output = argv[4];
    int nlr;
    FILE *fpw;
    FILE *fpm;

    // Allocate memory for input power spectrum and computation arrays
    double *wl=(double*) calloc(2*lmax+1,sizeof(double));
    double *lnjp1 = (double*) calloc(4 * lmax + 1, sizeof(double));
    double *lng = (double*) calloc(2 * lmax + 1, sizeof(double));
    double *g = (double*) calloc(2 * lmax + 1, sizeof(double));
    double *one_g = (double*) calloc(2 * lmax + 1, sizeof(double));
    double *one_jp1 = (double*) calloc(4 * lmax + 1, sizeof(double));
    double* m=(double*) calloc((lmax+1)*(lmax+1),sizeof(double));

    //Read in power spectrum 
    fpw=fopen(argv[1],"rb");
    nlr=fread((void*)wl,sizeof(double),2*lmax+1,fpw);
    if(nlr!=(2*lmax+1)){
        printf("Error reading power spectrum.\n");
        fclose(fpw);    
        return -1;  
    }   
    fclose(fpw);    

    // Precompute values for lnjp1, lng, g, and one_g
    for (int i = 0; i < 4 * lmax + 1; i++) {
        lnjp1[i] = log(i + 1.0);
        one_jp1[i] = 1.0 / (i + 1.0);
    }
    
    lng[0] = 0.0;
    for (int i = 1; i < 2 * lmax + 1; i++) {
        lng[i] = lng[i - 1] + log((i - 0.5) / i);
        one_g[i] = exp(-lng[i]);
    }

    for (int i = 0; i < 2 * lmax + 1; i++) {
        g[i] = exp(lng[i]);
    }

    double j1, j2, j3, j3_sum,pref;
    int J, J1minus, J2minus, J3minus;

    clock_t start, end;
    double cpu_time_used;

    // Calculating matrix M_l1l2 , dimensions j1xj2 = lmax * lmax
    start = clock();
    
    #pragma omp parallel for private(j1, j2, j3, J, J1minus, J2minus, J3minus, pref, j3_sum) \
    shared(m, wl, one_jp1, g, one_g) \
    schedule(dynamic)
    for(int i = 2; i < lmax + 1 ; i++) {
        for(int k = i;  k < lmax + 1; k++){
            j1 = (double)i;
            j2 = (double)k;
            j3_sum = 0.;
           
            for(int l = abs(i-k); l < i+k+1 ; l++){ 
                
                j3 = (double)l;
                J = i + k + l;
                if (J % 2 == 0) {
                    J1minus = -i + k + l;
                    J2minus = i - k + l;    
                    J3minus = i + k - l;
                    double threej_000 = calculate_threej_000(J, J1minus, J2minus, J3minus, one_jp1, g, one_g);
                    pref = ((2 * l) + 1) * wl[l];
                    j3_sum+= (pref * threej_000 * threej_000);
                }
            }
            m[i*(lmax+1)+k] =  j3_sum *  M_1_PI * 0.25;
      
        }
    }
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < lmax + 1; i++) {
        for (int j = i ; j < lmax + 1; j++) {  // j = i + 1 avoids unnecessary self-copy
            m[j * (lmax + 1) + i] = m[i * (lmax + 1) + j];  // Copy upper triangle to lower triangle
        }
    }
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < lmax + 1; i++) {
        for (int j = 0; j < lmax + 1; j++) {
            m[i * (lmax + 1) + j] *= (2. * j + 1.);
        }
    }

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("CPU time used : %g\n", cpu_time_used);
    
    // Write to file only if user requested
    if (strcmp(write_output, "yes") == 0) {
        FILE *fpm = fopen(argv[2], "wb");
        if (!fpm) {
            perror("Error opening output file");
            free(wl);
            free(m);
            return 1;
            }
        fwrite((void*)m, sizeof(double), (lmax+1)*(lmax+1), fpm);
        fclose(fpm);
        printf("Matrix written to %s.\n", argv[2]);
        }

    // Free allocated memory
    free(lnjp1);
    free(lng);
    free(g);
    free(one_g);
    free(one_jp1);

    return 0;
}

