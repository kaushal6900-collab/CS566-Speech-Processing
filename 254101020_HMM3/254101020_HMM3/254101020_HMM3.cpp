#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ---------- Read initial model ----------
void readInitialModel(const char *filename, double A[5][5], double B[5][32], double pi[5]) {
    FILE *f = fopen(filename,"r");
    char line[256];
	
    while(fgets(line, sizeof(line), f)) if (strstr(line, "A Matrix")) break;
    for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) fscanf(f, "%lf", &A[i][j]);

    while (fgets(line, sizeof(line), f)) if (strstr(line, "B Matrix")) break;
    for (int i = 0; i < 5; i++) for (int j = 0; j < 32; j++) fscanf(f, "%lf", &B[i][j]);

    while (fgets(line, sizeof(line), f)) if (strstr(line, "PI Matrix")) break;
    for (int i = 0; i < 5; i++) fscanf(f, "%lf", &pi[i]);

    fclose(f);
}

// ---------- Read observation sequence ----------
int readO(int O[200], const char *filename) {
    FILE *f = fopen(filename, "r");
    int T = 0;
    while (fscanf(f, "%d", &O[T]) == 1) {
        O[T] -= 1; // assumes input 1–32 → convert to 0–31
        T++;
    }
    fclose(f);
    return T;
}

// ---------- Forward algorithm ----------
void forward(int O[], int T, double A[5][5], double B[5][32], double pi[5],
             double alpha[200][5], double c[200]) {
    c[0] = 0.0;
    for (int i = 0; i < 5; i++) {
        alpha[0][i] = pi[i] * B[i][O[0]];
        c[0] += alpha[0][i];
    }
    if (c[0] == 0) c[0] = 1e-30;
    c[0] = 1.0 / c[0];
    for (int i = 0; i < 5; i++) alpha[0][i] *= c[0];

    for (int t = 1; t < T; t++) {
        c[t] = 0.0;
        for (int i = 0; i < 5; i++) {
            alpha[t][i] = 0.0;
            for (int j = 0; j < 5; j++) alpha[t][i] += alpha[t - 1][j] * A[j][i];
            alpha[t][i] *= B[i][O[t]];
            c[t] += alpha[t][i];
        }
        if (c[t] == 0) c[t] = 1e-30;
        c[t] = 1.0 / c[t];
        for (int i = 0; i < 5; i++) alpha[t][i] *= c[t];
    }
}

// ---------- Backward algorithm ----------
void backward(int O[], int T, double A[5][5], double B[5][32], double beta[200][5], double c[200]) {
    for (int i = 0; i < 5; i++) beta[T - 1][i] = c[T - 1];
    for (int t = T - 2; t >= 0; t--) {
        for (int i = 0; i < 5; i++) {
            beta[t][i] = 0.0;
            for (int j = 0; j < 5; j++) beta[t][i] += A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
            beta[t][i] *= c[t];
        }
    }
}

// ---------- Gamma and Xi ----------
void computeGammaXi(int O[], int T, double A[5][5], double B[5][32],
                    double alpha[200][5], double beta[200][5],
                    double gamma[200][5], double xi[200][5][5]) {
    for (int t = 0; t < T - 1; t++) {
        double denom = 0.0;
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                denom += alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
        if (denom == 0) denom = 1e-30;

        for (int i = 0; i < 5; i++) {
            gamma[t][i] = 0.0;
            for (int j = 0; j < 5; j++) {
                xi[t][i][j] = (alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j]) / denom;
                gamma[t][i] += xi[t][i][j];
            }
        }
    }
    for (int i = 0; i < 5; i++) gamma[T - 1][i] = alpha[T - 1][i];
}

// ---------- Normalize A ----------
void renormA(double A[5][5]) {
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++)
            if (!(j == i || j == i + 1)) A[i][j] = 0.0;
        if (i == 4) { for (int j = 0; j < 4; j++) A[4][j] = 0.0; A[4][4] = 1.0; continue; }

        double sum = A[i][i] + A[i][i + 1];
        if (sum == 0) { A[i][i] = 0.9; A[i][i + 1] = 0.1; sum = 1.0; }
        A[i][i] /= sum; A[i][i + 1] /= sum;
    }
}

// ---------- Normalize B (professor's exact order) ----------
void renormB(double B[5][32]) {
    const double floorv = 1e-30;
    for (int i = 0; i < 5; i++) {
        double sum = 0.0;
        for (int k = 0; k < 32; k++) sum += B[i][k];
        if (sum == 0) sum = 1e-30;
        for (int k = 0; k < 32; k++) B[i][k] /= sum;

        for (int k = 0; k < 32; k++) if (B[i][k] < floorv) B[i][k] = floorv;

        double sum2 = 0.0;
        for (int k = 0; k < 32; k++) sum2 += B[i][k];
        for (int k = 0; k < 32; k++) B[i][k] /= sum2;
    }
}

// ---------- Re-estimation ----------
void reestimate(int O[], int T, double A[5][5], double B[5][32],
                double pi[5], double gamma[200][5], double xi[200][5][5]) {
    for (int i = 0; i < 5; i++) pi[i] = gamma[0][i];

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            double num = 0, den = 0;
            for (int t = 0; t < T - 1; t++) { num += xi[t][i][j]; den += gamma[t][i]; }
            if (den == 0) den = 1e-30;
            A[i][j] = num / den;
        }
    }

    for (int i = 0; i < 5; i++) {
        for (int k = 0; k < 32; k++) {
            double num = 0, den = 0;
            for (int t = 0; t < T; t++) { if (O[t] == k) num += gamma[t][i]; den += gamma[t][i]; }
            if (den == 0) den = 1e-30;
            B[i][k] = num / den;
        }
    }

    renormA(A);
    renormB(B);
}

// ---------- Log-likelihood ----------
double computeLogLikelihood(double c[200], int T) {
    double logP = 0.0;
    for (int t = 0; t < T; t++) logP += log(1.0 / c[t]);
    return logP;
}

// ---------- Training ----------
void trainHMM(int O[], int T, double A[5][5], double B[5][32], double pi[5]) {
    double alpha[200][5], beta[200][5], gamma[200][5], xi[200][5][5];
    double c[200];
    int maxIter = 20;

    for (int iter = 0; iter < maxIter; iter++) {
        forward(O, T, A, B, pi, alpha, c);
        backward(O, T, A, B, beta, c);
        computeGammaXi(O, T, A, B, alpha, beta, gamma, xi);
        reestimate(O, T, A, B, pi, gamma, xi);
        printf("Iteration %2d : log[P(O|λ)] = %lf\n", iter + 1, computeLogLikelihood(c, T));
    }
}

// ---------- Save model ----------
void saveModel(const char *aFile, const char *bFile, const char *piFile,
               double A[5][5], double B[5][32], double pi[5]) {
    FILE *fa = fopen(aFile, "w");
    for (int i = 0; i < 5; i++) { for (int j = 0; j < 5; j++) fprintf(fa, "%.6le\t", A[i][j]); fprintf(fa, "\n"); }
    fclose(fa);

    FILE *fb = fopen(bFile, "w");
    for (int i = 0; i < 5; i++) { for (int j = 0; j < 32; j++) fprintf(fb, "%.6le\t", B[i][j]); fprintf(fb, "\n"); }
    fclose(fb);

    FILE *fp = fopen(piFile, "w");
    for (int i = 0; i < 5; i++) fprintf(fp, "%.6le\t", pi[i]);
    fprintf(fp, "\n");
    fclose(fp);
}

// ---------- Main ----------
int _tmain(int argc, _TCHAR* argv[]) {
    double A[5][5], B[5][32], pi[5];
    int O[200];
    readInitialModel("Initial_Model.txt", A, B, pi);
    int T = readO(O, "HMM_OBSERVATION_SEQUENCE_1.txt");
    trainHMM(O, T, A, B, pi);
    saveModel("HMM_AIJ_FINAL.txt", "HMM_BJK_FINAL.txt", "HMM_PI_FINAL.txt", A, B, pi);
    return 0;
}