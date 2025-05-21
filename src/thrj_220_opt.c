#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

int lmax;


double calculate_threej_000(int J, int J1minus, int J2minus, int J3minus, double* one_jp1, double* g, double* one_g) {
    if (J % 2 == 0) {
        //return pow(-1, J / 2) * sqrt(one_jp1[J] * g[J1minus / 2] * g[J2minus / 2] * g[J3minus / 2] * one_g[J / 2]);
        int half_J1minus = J1minus / 2;
        int half_J2minus = J2minus / 2;
        int half_J3minus = J3minus / 2;
        int half_J = J / 2;
        //return sqrt(one_jp1[J] * g[J1minus / 2] * g[J2minus / 2] * g[J3minus / 2] * one_g[J / 2]);
        return pow(-1, J / 2) * sqrt(one_jp1[J] * g[half_J1minus] * g[half_J2minus] * g[half_J3minus] * one_g[half_J]);
    }else{
        printf("Error: J is odd\n");
        return 0.0;
    }
}

// Function to calculate 3j (0,-2,2) symbol depending on even or odd J
double calculate_threej_022(double j1, double j2, double j3, double* one_jp1, double* g, double* one_g) {
    double threej_022 = 0.0;
    double J = j1 + j2 + j3;
    double J1minus = -j1 + j2 + j3;
    double J2minus = j1 - j2 + j3;
    double J3minus = j1 + j2 - j3;
    if ((int)J % 2 == 0) {  // J is even


        double pref_1 =  1. / sqrt((j2-1.)*(j2+2.)*(j3-1.)*j3);
        double pref_2_1 = sqrt(j2*(j2+1.)*(j3+1.)*(j3+2.));
        double pref_2 = 2.*sqrt((j2+1)*(j3+1)*(j3+2)/j2);
        double pref_3 = 1. - (0.5*((J+2)*(J1minus+1.)/((j2+1.)*(j3+1.))));
        double pref_5 = 0.5 * sqrt((J+2.)*(J1minus+1.)*(J2minus+1.)*(J3minus)*(J+3.)*(J1minus+2.)*(J2minus+2.)*(J3minus-1.)/(j2*(j2+1.)*(j3+1.)*(j3+2.)));

        double threej_000 = calculate_threej_000((int)J, (int)J1minus, (int)J2minus, (int)J3minus, one_jp1, g, one_g);
        
        double threej_000_2 = calculate_threej_000((int)J+2, (int)J1minus+2, (int)J2minus+2, (int)J3minus-2, one_jp1, g, one_g);

        threej_022 = pref_1 * (((pref_2_1 + (pref_2 * pref_3))*threej_000) + (pref_5 * threej_000_2));
    } else {  // J is odd

        double pref_0 = sqrt((J + 2.0) * (J1minus + 1.0) * (J2minus + 1.0) * J3minus) / sqrt((j2 - 1.0) * (j2 + 2.0) * (j3 - 1.0) * j3);
        double pref_1 = sqrt((j3 + 2.0) / (j2 * (j2 + 1.0) * (j3 + 1.0)));
        double pref_2 = sqrt((j2 + 1.0) * (j3 + 2.0) / (j2 * (j3 + 1.0)));
        double pref_3 = -0.5 * sqrt((J + 3.0) * (J + 4.0) * (J1minus + 2.0) * (J1minus + 3.0) / (j2 * (j2 + 1.0) * (j3 + 1.0) * (j3 + 2.0)));

        //double threej_000_1 = calculate_threej_000(j1, j2, j3 + 1, one_jp1, g, one_g);
        double threej_000_1 = calculate_threej_000((int)J+1, (int)J1minus+1, (int)J2minus+1, (int)J3minus-1, one_jp1, g, one_g);
        //double threej_000_2 = calculate_threej_000(j1, j2 + 1, j3 + 2, one_jp1, g, one_g);
        double threej_000_2 = calculate_threej_000((int)J+3, (int)J1minus+3, (int)J2minus+1, (int)J3minus-1, one_jp1, g, one_g);

        threej_022 = pref_0 * ((-1.0 * (pref_1 + pref_2) * threej_000_1) + (pref_3 * threej_000_2));
        //threej_022 = -1.0 * pref_0 * (( -1.0 * (pref_1 + pref_2) * threej_000_1) + (pref_3 * threej_000_2));
    }
    return threej_022;
    //return (pow(-1, - j1 + j2 - 2) * threej_022);
}

int main(int argc, char *argv[]) {
     
    lmax=atoi(argv[3]);
    int nlr;
    char *pol = argv[4];
    FILE *fpw;
    FILE *fpm;

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
    
    lng[0] = 0.0;
    for (int i = 1; i < 2 * lmax + 1; i++) {
        lng[i] = lng[i - 1] + log((i - 0.5) / i);
        //g[i] = exp(lng[i]);
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
    // For j1 > j2
    // Start at 2 as -l1 <=  m1 <= l1
    //#pragma omp parallel for schedule(dynamic) private(j1, j2, j3, j3_sum, pref)
    for(int i = 2; i < lmax + 1 ; i++) {
        for(int k = i;  k < lmax + 1; k++){
            j1 = (double)i;
            j2 = (double)k;

            j3_sum = 0.;

            for(int l = abs(i-k); l < i+k+1 ; l++){  //is this right for j3????
                j3 = (double)l;
                //printf("ORIGINAL ---> j1, j2, j3: %d, %d, %d\n", i, k, l);

                //Permutation
                int tempj1 = i; // Store j1 temporarily
                j1 = l;   
                j3 = k;   
                j2 = tempj1; 

                if(j3>=2){
                    
                    double threej_022 = calculate_threej_022(j1, j2, j3, one_jp1, g, one_g);
        
                    // positive sign for KEE, negative sign for KEB
                    //if(strcmp(pol, "EE") == 0){
                    //int parity = (int)(j1 + j2 + j3) % 2;
                    //pref = ((2 * l) + 1) * wl[l] * (parity == 0 ? 4 : 0);
                    pref = ((2 * l) + 1) * wl[(int)l] * pow((1 + pow(-1,(j1+j2+j3))),2);
                    j3_sum+= (pref * threej_022 * threej_022);
                }
            }
            
            m[i*(lmax+1)+k] =  j3_sum *  M_1_PI / 16.0;
        }
    }

    for(int i=0 ; i<lmax+1 ; i++) {
        for (int j=i ; j < lmax+1 ; j++) {
            m[j*(lmax+1)+i]=m[i*(lmax+1)+j];
        }
    }
    for(int i=0;i<lmax+1;i++) {
        for (int j=0;j<lmax+1;j++) {
            m[i*(lmax+1)+j]*=(2.*j+1.);  // is it quicker to multiply the entire row as they all have the same l2 ? 
        }
    }

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("CPU time used : %g\n", cpu_time_used);

    fpm=fopen(argv[2],"wb");
    if (fpm == NULL) {
        perror("Error opening file");
        return 1;
    }
    fwrite((void*)m,sizeof(double),(lmax+1)*(lmax+1),fpm);
    fclose(fpm);

    // Free allocated memory
    free(lnjp1);
    free(lng);
    free(g);
    free(one_g);
    free(one_jp1);

    return 0;
}

