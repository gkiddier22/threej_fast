/*
 * thrj_220_opt.c
 * 
 * Description: Optimized C code for calculating the coupling matrix for CMB polarisation maps
 * Author: Georgia Kiddier
 * Date: 18-02-2025
 * License: MIT (or specify your chosen license)
 * 
 * Compilation:
 *     /opt/homebrew/bin/gcc-14 -O3 -march=native -o thrj_220_opt thrj_220_opt.c -lm -ffast-math
 * 
 * Usage:
 *     ./thrj_220_opt <power_spectrum_file> <output_matrix_file> <lmax>
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

/**
 * @brief Computes the Wigner 3j symbol for the configuration (0, 0, 0).
 *
 * This function evaluates the Wigner 3j symbol:
 * \[
 * \begin{pmatrix}
 * j_1 & j_2 & j_3 \\
 * 0   & 0   & 0
 * \end{pmatrix}
 * \]
 * using precomputed arrays to minimize computational cost. It assumes that the triangle
 * conditions are satisfied and that the total angular momentum quantum number \( J = j_1 + j_2 + j_3 \)
 * is even, which is a necessary condition for the symbol to be nonzero.
 *
 * @param J         The sum \( j_1 + j_2 + j_3 \), assumed to be even.
 * @param J1minus   The value \( -j_1 + j_2 + j_3 \).
 * @param J2minus   The value \( j_1 - j_2 + j_3 \).
 * @param J3minus   The value \( j_1 + j_2 - j_3 \).
 * @param one_jp1   Array containing \( \frac{1}{2j+1} \) values for normalization.
 * @param g         Array containing factorial-related values: \( g(n) = n! \).
 * @param one_g     Array containing \( \frac{1}{g(n)} = \frac{1}{n!} \).
 *
 * @return The computed Wigner 3j symbol for (0, 0, 0).
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


/**
 * @brief Computes the Wigner 3j symbol for the configuration (0, -2, 2)
 *        using a fast approximation based on precomputed arrays and recursion.
 *
 * This function approximates the Wigner 3j symbol:
 * \[
 * \begin{pmatrix}
 * j_1 & j_2 & j_3 \\
 * 0   & -2  & 2
 * \end{pmatrix}
 * \]
 * by relating it to the analytically simpler (0, 0, 0) configuration, with corrections
 * based on recursion relations and precomputed normalization arrays. It relies on 
 * efficient reuse of cached factorial terms and normalization factors for performance.
 *
 * @param j1        Angular momentum quantum number \( j_1 \)
 * @param j2        Angular momentum quantum number \( j_2 \)
 * @param j3        Angular momentum quantum number \( j_3 \)
 * @param J         The total \( J = j_1 + j_2 + j_3 \)
 * @param J1minus   \( -j_1 + j_2 + j_3 \)
 * @param J2minus   \( j_1 - j_2 + j_3 \)
 * @param J3minus   \( j_1 + j_2 - j_3 \)
 * @param one_jp1   Array of \( \frac{1}{2j + 1} \) for normalization
 * @param g         Array of factorial values: \( g(n) = n! \)
 * @param one_g     Array of inverse factorials: \( \frac{1}{n!} \)
 * @param one_r_j   Array of \( \frac{1}{\sqrt{J}} \), used for normalization in recursion terms
 *
 * @return Approximated Wigner 3j symbol for (0, -2, 2)
 */
double calculate_threej_022_alt(int j1, int j2, int j3, int J, int J1minus, int J2minus, int J3minus, double* one_jp1, double* g, double* one_g, double* one_r_j) {
    double threej_022 = 0.0;
    double lmbda = sqrt(j2*(j2+1.)*(j3+1.)*(j3+2.));


    int J2MP1 = J2minus + 1;
    int J3MM1 = J3minus - 1;

        double pref_1 = one_r_j[(j2-1)] *one_r_j[(j2+2)] * one_r_j[(j3-1)] * one_r_j[j3];
        
        double pref_2_1 = lmbda;
        
        double pref_2 = 2. * lmbda * one_jp1[(j2-1)];
   
        double lmbda2 = (J+2)*(J1minus+1);
   
        double pref_3 = 1. - 0.5 * lmbda2 * one_jp1[j2] * one_jp1[j3];
        
        double pref_5 = 0.5 * sqrt(lmbda2*(J2MP1)*(J3minus)*(J+3.)*(J1minus+2.)*(J2minus+2.)*(J3MM1)) / lmbda;
       
        double threej_000 = calculate_threej_000(J, J1minus, J2minus, J3minus, one_jp1, g, one_g);
        
        double threej_000_2 = calculate_threej_000(J+2, J1minus+2, J2minus+2, J3minus-2, one_jp1, g, one_g);

        threej_022 = pref_1 * (((pref_2_1 + (pref_2 * pref_3))*threej_000) + (pref_5 * threej_000_2));

    return threej_022;
    
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <input_file> <output_file> <lmax> <pol> <write_output (yes/no)>\n", argv[0]);
        return 1;
    }
     
    lmax=atoi(argv[3]);
    int nlr;
    char *pol = argv[4];
    char *write_output = argv[5];
    FILE *fpw;
    FILE *fpm;

    double *wl=(double*) calloc(2*lmax+1,sizeof(double));
    double *lnjp1 = (double*) calloc(4 * lmax + 1, sizeof(double));
    double *lng = (double*) calloc(2 * lmax + 1, sizeof(double));
    double *g = (double*) calloc(2 * lmax + 1, sizeof(double));
    double *one_g = (double*) calloc(2 * lmax + 1, sizeof(double));
    double *one_jp1 = (double*) calloc(4 * lmax + 1, sizeof(double));
    double *one_r_j = (double*) calloc((2*lmax)*(2*lmax), sizeof(double));
    double* m=(double*) calloc((lmax+1)*(lmax+1),sizeof(double));

    if (!wl || !lnjp1 || !lng || !g || !one_g || !one_jp1 || !one_r_j || !m) {
    fprintf(stderr, "Memory allocation failed!\n");
    return -1;
    }

    //Read in power spectrum 
    fpw=fopen(argv[1],"rb");
    nlr=fread((void*)wl,sizeof(double),2*lmax+1,fpw);
    if(nlr!=(2*lmax+1)){
        printf("Error reading power spectrum.\n");
        fclose(fpw);    
        return -1;  
    }   
    fclose(fpw);    
    //printf("lmax*lmax: %d\n", lmax*lmax);

    // Validate the mode
    if (strcmp(pol, "EE") != 0 && strcmp(pol, "EB") != 0) {
        fprintf(stderr, "Error: pol must be either 'EE' or 'EB'.\n");
        return 1;
    }

    // Precompute values for lnjp1, lng, g, and one_g
    for (int i = 0; i < 4 * lmax + 1; i++) {
        lnjp1[i] = log(i + 1.0);
        one_jp1[i] = 1.0 / (i + 1.0);
    }

    one_r_j[0] = 0.0;
    for (int i = 1; i < lmax*lmax; i++) {
        //cache 1/sqrt(J)
        one_r_j[i] = 1.0 / sqrt(i);
    }
    
    lng[0] = 0.0;
    for (int i = 1; i < 2 * lmax + 1; i++) {
        lng[i] = lng[i - 1] + log((i - 0.5) / i);
        one_g[i] = exp(-lng[i]);
    }

    for (int i = 0; i < 2 * lmax + 1; i++) {
        g[i] = exp(lng[i]);
    }

    #define MAX(a, b) ((a) > (b) ? (a) : (b))

    double j1, j2, j3, j3_sum,pref;
    int J, J1minus, J2minus, J3minus;
    int size = (lmax + 1) * (lmax + 1);

    clock_t start, end;
    double cpu_time_used;

    // Calculating matrix M_l1l2 , dimensions j1xj2 = lmax * lmax
    
    start = clock();
    double scale_factor = M_1_PI * 0.0625;

    #pragma omp parallel for private(j1, j2, j3, J, J1minus, J2minus, J3minus, pref, j3_sum) \
    shared(m, wl, one_jp1, g, one_g) \
    schedule(dynamic)
    for(int i = 2; i < lmax + 1 ; i++) {
        for(int k = MAX(2, i);  k < lmax + 1; k++){

            j3_sum = 0.;
            int jbase = i + k;

            for(int l = abs(i-k); l < jbase+1 ; l++){  
                // j3 = (double)l;
                double threej_022 = 0.0;
                //Permutation necessary as we need (-2,2,0) symbol; function computes (0,-2,2)

                    J = jbase + l;
                    //Only get symbol if J is even
                    if (J % 2 == 0) {
                        J1minus = -l + jbase;
                        J2minus = l - i + k;
                        J3minus = l + i - k;
                    
                        threej_022 = calculate_threej_022_alt(l, i, k, J, J1minus, J2minus, J3minus, one_jp1, g, one_g, one_r_j);
                    
                        pref = ((2 * l) + 1) * wl[l] * 4;
                        j3_sum+= (pref * threej_022 * threej_022);
                    
                    }
            }
            m[i*(lmax+1)+k] =  j3_sum * scale_factor;
      
        }
    }


    #pragma omp parallel for collapse(2)
    for (int i = 0; i < lmax + 1; i++) {
        for (int j = i + 1; j < lmax + 1; j++) {  // j = i + 1 avoids unnecessary self-copy
            m[j * (lmax + 1) + i] = m[i * (lmax + 1) + j];  // Copy upper triangle to lower triangle
        }
    }
    #pragma omp parallel for collapse(2)
    for (int idx = 0; idx < size; idx++) {
        int j = idx % (lmax + 1);
        m[idx] *= (2. * j + 1.);
    }


    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;


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
    free(one_r_j);

    return 0;
}