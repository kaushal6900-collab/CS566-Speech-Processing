#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <windows.h>

#define MAX_SAMPLES 25000
#define NORM_RANGE 100
#define FRAME_SIZE 320
#define STEADY_START 5000
#define NUM_FRAMES 5
#define P 12
#define Q 12
#define PI (22.0/7.0)

// Pre-defined weights for the Tokhura distance calculation.
double Tokhuraweights[Q] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

//  DC offset removal
void dcshift(double *samples, int n) {
    double sum = 0.0;
    for(int i=0; i<n; i++) {
        sum = sum + samples[i];
    }
    double mean = sum/n;

    for(int i=0; i<n; i++) {
        samples[i] = samples[i] - mean; 
    }
}

// Normalization process
void normalization(double *sample, int n) {
    double max = 0.0;

    for(int i=0; i<n; i++) {
        if (fabs(sample[i]) > max) {
            max = fabs(sample[i]);
        }
    }

    for(int i=0; i<n; i++) {
        sample[i] = NORM_RANGE * sample[i] / max;
    }
}

// Extracting steady part of the audio signal.
void steadyframes(double *sample, double frames[NUM_FRAMES][FRAME_SIZE]) {
    for(int i=0; i<NUM_FRAMES; i++) {
        int offset = STEADY_START + i * FRAME_SIZE;
        for(int j=0; j<FRAME_SIZE; j++) {
            frames[i][j] = sample[offset + j];
        }
    }
}

//  Hamming window implementation
void hammingwindow(double frame[]) {
    for(int i=0; i<FRAME_SIZE; i++) {
        double w = 0.54 - 0.46 * cos(2 * PI * i / (FRAME_SIZE - 1));
        frame[i] *= w;
    }
}

// Calculation of autocorrelation coefficients 
void autocorrelation(double frame[], double R[]) {
     for(int i=0; i<=P; i++) {
         R[i] = 0.0;
         for(int j=0; j<FRAME_SIZE-i; j++) {
             R[i] += frame[j] * frame[j+i];
         }
    }
}

// Implementation of Levinson-Durbin algorithm 
void levinson_durbin_Algo(double R[], double A[]) {
    double E = R[0];
    double a[P] = {0.0};

    for(int m=1; m<=P; m++) {
        double k = -R[m];
        for(int i=0; i<m-1; i++) {
            k -= a[i] * R[m-i-1];
        }
        k = k / E;

        double a_prev[P] = {0.0};
        for(int i=0; i<m-1; i++) {
            a_prev[i] = a[i];
        }

        for(int i=m-2; i>=0; i--) {
            a[i] = a_prev[i] + k * a_prev[m-2-i];
        }
        a[m-1] = k;

        E = E * (1 - k*k);
    }

    for(int i=0; i<P; i++) {
        A[i] = -a[i];
    }
}

// Converts the LPC coefficients (A) into Cepstral coefficients (C).
void cepstral(double A[], double R0, double C[]) {
    C[0] = log(R0 * R0);

    for(int m=1; m<=Q; m++) {
        double sum = 0.0;
        for(int k=1; k<m; k++) {
            sum += ((double)k / m) * C[k] * A[m-k-1];
        }
        if (m <= P) { C[m] = A[m-1] + sum; }
        else { C[m] = sum; }
    }
}

// Calculation of weighted Tokhura distance 
double tokhuraDistance(double C1[Q+1], double C2[Q+1]) {
    double distance = 0.0;
    for(int m=1; m<=Q; m++) {
        double diff = C1[m] - C2[m];
        distance += Tokhuraweights[m-1] * diff * diff;
    }
    return distance;
}


int _tmain(int argc, _TCHAR* argv[])
{
    // Training start
    
    printf("training...");
    char inputfile[100];
    char outputfile[100];
    
    int num = 5;
    while(num--) {
        printf("...");
        char vowel = ' ';
        if (num == 4) {vowel = 'a';}
        if (num == 3) {vowel = 'e';}
        if (num == 2) {vowel = 'i';}
        if (num == 1) {vowel = 'o';}
        if (num == 0) {vowel = 'u';}
        
        double Ci_sum[NUM_FRAMES][Q] = {0};

        for(int x=1; x<21; x++) {
            double samples[MAX_SAMPLES];
    
            sprintf(inputfile, "vowels_train/%c/254101020_%c_%d.txt", vowel, vowel, x);
            FILE *f1 = fopen(inputfile, "r");
    
            char line[256];
            for(int y=0; y<5; y++) {
                fgets(line, sizeof(line), f1);
            }
    
            int n = 0;
            while(fscanf(f1, "%lf", &samples[n]) == 1) {
                n++;
            }

            dcshift(samples, n);
            normalization(samples, n);

            double frames[NUM_FRAMES][FRAME_SIZE];
            steadyframes(samples, frames);

            double R[P+1] = {0};
            double Rnorm[P+1];
            double A[P];
            double C[Q+1];

            for(int i=0; i<NUM_FRAMES; i++) {
                hammingwindow(frames[i]);
                autocorrelation(frames[i], R);

                for(int j=0; j<=P; j++) {
                    Rnorm[j] = R[j] / R[0];
                }

                levinson_durbin_Algo(Rnorm, A);
                cepstral(A, R[0], C);
                
                for(int m=1; m<=Q; m++) {
                    double w = 1 + (Q/2.0) * sin(PI * m / Q);
                    C[m] = C[m] * w;

                    Ci_sum[i][m-1] += C[m];
                }

            }

            fclose(f1);
        }
        

        sprintf(outputfile, "CiAvg/%c_CiAvg.txt", vowel);
        FILE *f2 = fopen(outputfile, "w");
        for(int i=0; i<NUM_FRAMES; i++) {
            for(int j=0; j<Q; j++) {
                double avg = Ci_sum[i][j] / 20.0;
                fprintf(f2, "%lf\n", avg);
            }
        }

        fclose(f2);
    }
    printf("  training done!\n");

    Sleep(2000);

    //  Testing start
    
    printf("\ntesting start\n\n");

    num = 5;
    while(num--) {
        char vowel = ' ';
        if (num == 4) {vowel = 'a';}
        if (num == 3) {vowel = 'e';}
        if (num == 2) {vowel = 'i';}
        if (num == 1) {vowel = 'o';}
        if (num == 0) {vowel = 'u';}

        for(int x=21; x<31; x++) {
            double samples[MAX_SAMPLES];
    
            sprintf(inputfile, "vowels_test/%c/254101020_%c_%d.txt", vowel, vowel, x);
            FILE *f1 = fopen(inputfile, "r");
    
            char line[256];
            for(int y=0; y<5; y++) {
                fgets(line, sizeof(line), f1);
            }
    
            int n = 0;
            while(fscanf(f1, "%lf", &samples[n]) == 1) {
                n++;
            }

            dcshift(samples, n);
            normalization(samples, n);

            double frames[NUM_FRAMES][FRAME_SIZE];
            steadyframes(samples, frames);

            double R[P+1] = {0};
            double Rnorm[P+1];
            double A[P];
            double C[Q+1];
            double C_test[NUM_FRAMES][Q+1];

            for(int i=0; i<NUM_FRAMES; i++) {
                hammingwindow(frames[i]);
                autocorrelation(frames[i], R);

                for(int j=0; j<=P; j++) {
                    Rnorm[j] = R[j] / R[0];
                }

                levinson_durbin_Algo(Rnorm, A);
                cepstral(A, R[0], C);
                
                for(int m=1; m<=Q; m++) {
                    double w = 1 + (Q/2.0) * sin(PI * m / Q);
                    C[m] = C[m] * w;
                }

                for(int m=1; m<=Q; m++) {
                    C_test[i][m] = C[m];
                }
            }

            char inputfile2[100];
            
            double min_distance = 1e9;
            char vowel_recognized = ' ';

            for(int y=0; y<5; y++) {
                char ref_vowel = ' ';
                if (y == 0) ref_vowel = 'a';
                if (y == 1) ref_vowel = 'e';
                if (y == 2) ref_vowel = 'i';
                if (y == 3) ref_vowel = 'o';
                if (y == 4) ref_vowel = 'u';

                sprintf(inputfile2, "CiAvg/%c_CiAvg.txt", ref_vowel);
                FILE *f2 = fopen(inputfile2, "r");

                double C_ref[NUM_FRAMES][Q+1];
                for(int i=0; i<NUM_FRAMES; i++) {
                    for(int j=1; j<=Q; j++) {
                        fscanf(f2, "%lf", &C_ref[i][j]);
                    }
                }

                double total_distance = 0.0;
                for(int i=0; i<NUM_FRAMES; i++) {
                    total_distance += tokhuraDistance(C_test[i], C_ref[i]);
                }
                double avg_distance = total_distance / NUM_FRAMES;

                if (avg_distance < min_distance) {
                    min_distance = avg_distance;
                    vowel_recognized = ref_vowel;
                }

                fclose(f2);
            
            }
            printf("File %c_%d recognized as: %c\n", vowel, x, vowel_recognized);

            fclose(f1);
        }
        printf("\n");
    }

    return 0;
}