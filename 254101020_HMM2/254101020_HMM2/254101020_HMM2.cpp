#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Reads final trained HMM matrices A (5x5) and B (5x32)
void readFinalAB(const char *aFile, const char *bFile, double A[5][5], double B[5][32]) {
    FILE *fa = fopen(aFile,"r");
    FILE *fb = fopen(bFile,"r");

    // Read A matrix (5x5)
    for(int i=0; i<5; i++) {
        for(int j=0; j<5; j++) {
            fscanf(fa, "%lf", &A[i][j]);
        }
    }

    // Read B matrix (5x32)
    for(int i=0; i<5; i++) {
        for(int j=0; j<32; j++) {
            fscanf(fb, "%lf", &B[i][j]);
        }
    }

	fclose(fa);
    fclose(fb);
}


void readPi(const char *filename, double pi[5]) {
    FILE *f = fopen(filename,"r");
    char line[1024];

    // Find the line that contains "PI Matrix"
    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, "PI Matrix")) {
            break;
        }
    }

    // Read 5 values of π
    for (int i = 0; i < 5; i++) {
        fscanf(f, "%lf", &pi[i]);
    }

    fclose(f);
}


// reads the observation sequence from file into array O[]
int readO(int O[200], const char *filename) {
	FILE *f = fopen(filename, "r");
	int T = 0;
	while(fscanf(f, "%d", &O[T]) == 1) {
		T++;
	}
	fclose(f);
	return T;
}

// Viterbi algorithm for HMM with 5 states and up to 200 observations
double viterbiAlgorithm(double A[5][5], double B[5][32], double pi[5], int O[200], int T) {
    double delta[200][5];
    int psi[200][5];
    int Q[200];

    // 1. Initialization
    for(int i=0; i<5; i++) {
        delta[0][i] = pi[i] * B[i][O[0]-1];
        psi[0][i] = 0;
    }

    // 2. Recursion
    for(int t=1; t<T; t++) {
        for(int j=0; j<5; j++) {
            double maxVal = 0.0;
            int maxState = 0;
            for(int i=0; i<5; i++) {
                double val = delta[t-1][i] * A[i][j];
                if (val > maxVal) {
                    maxVal = val;
                    maxState = i;
                }
            }
            delta[t][j] = maxVal * B[j][O[t]-1];
            psi[t][j] = maxState;
        }
    }

    // 3. Termination
    double bestProb = 0.0;
    int bestLastState = 0;
    for(int i=0; i<5; i++) {
        if (delta[T-1][i] > bestProb) {
            bestProb = delta[T-1][i];
            bestLastState = i;
        }
    }

    // 4. Backtracking
    Q[T-1] = bestLastState;
    for(int t=T-2; t>=0; t--) {
        Q[t] = psi[t+1][Q[t+1]];
    }

    // 5. Output
    printf("Most probable state sequence:\n");
    for(int t=0; t<T; t++)
        printf("%d ", Q[t]+1);  // +1 for 1-based state indexing
    
	return bestProb;
}


int _tmain(int argc, _TCHAR* argv[])
{
	double A[5][5], B[5][32], pi[5];
	int O[200];
	
	readFinalAB("HMM_AIJ_FINAL.txt", "HMM_BJK_FINAL.txt", A, B);
    readPi("Initial_Model.txt", pi);
	int T = readO(O, "HMM_OBSERVATION_SEQUENCE_1.txt");

	double bestProb = viterbiAlgorithm(A, B, pi, O, T);
	printf("\nBest path probability = %.10e\n", bestProb);

	return 0;
}