#include "stdafx.h" 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define MAX_SAMPLES 25000
#define FRAME_SIZE 320
#define LPC_ORDER 12
#define NUM_STEADY_FRAMES 5
#define NUM_VOWELS 5
#define NUM_UTTERANCES 35

const char* rollNumber = "254101020"; 

void processFile(const char* inputFilename, const char* outputFilename);
void removeDCOffset(double* signal, int numSamples);
void normalizeSignal(double* signal, int numSamples);
void calculateAutocorrelation(const double* frame, double* R);
void levinsonDurbin(const double* R, double* A, double* E);
int findSteadyFrames(const double* signal, int numSamples, int* startFrameIndex);

// Loops through all sound files and processes them.
int main() {
    char vowels[NUM_VOWELS] = {'a', 'e', 'i', 'o', 'u'};
    char inputFilename[100];
    char outputFilename[100];
    int filesProcessed = 0;

    for (int i = 0; i < NUM_VOWELS; i++) {
        for (int j = 1; j <= NUM_UTTERANCES; j++) {
            sprintf(inputFilename, "%s_%c_%d.txt", rollNumber, vowels[i], j);
            sprintf(outputFilename, "%s_%c_%d_ai.txt", rollNumber, vowels[i], j);
            
            processFile(inputFilename, outputFilename);
            filesProcessed++;
        }
    }
    
    printf("\nFinished. The program attempted to process a total of %d files.\n", filesProcessed);
    
    printf("Press Enter to exit...");
    getchar(); 
    return 0;
}

// Processes one audio file to find its LPC coefficients.
void processFile(const char* inputFilename, const char* outputFilename) {
    FILE *inputFile = fopen(inputFilename, "r");
    if (!inputFile) {
        fprintf(stderr, "Warning: Could not open %s. Skipping.\n", inputFilename);
        return;
    }

    double signal[MAX_SAMPLES];
    int numSamples = 0;
    while (numSamples < MAX_SAMPLES && fscanf(inputFile, "%lf", &signal[numSamples]) == 1) {
        numSamples++;
    }
    fclose(inputFile);

    if (numSamples < FRAME_SIZE * NUM_STEADY_FRAMES) {
        fprintf(stderr, "Warning: Not enough samples in %s. Skipping.\n", inputFilename);
        return;
    }
    
    removeDCOffset(signal, numSamples);
    normalizeSignal(signal, numSamples);

    int steadyPartStartIndex;
    if (!findSteadyFrames(signal, numSamples, &steadyPartStartIndex)) {
        fprintf(stderr, "Warning: Could not find a stable region in %s. Skipping.\n", inputFilename);
        return;
    }

    FILE *outputFile = fopen(outputFilename, "w");
    if (!outputFile) {
        fprintf(stderr, "Error: Could not create output file %s.\n", outputFilename);
        return;
    }

    for (int i = 0; i < NUM_STEADY_FRAMES; i++) {
        int frameStartIndex = steadyPartStartIndex + (i * FRAME_SIZE);
        double *currentFrame = &signal[frameStartIndex];
        double R[LPC_ORDER + 1], A[LPC_ORDER + 1], E;

        calculateAutocorrelation(currentFrame, R);
        levinsonDurbin(R, A, &E);

        for (int j = 1; j <= LPC_ORDER; j++) {
            fprintf(outputFile, "%lf\n ", A[j]);
        }
        fprintf(outputFile, "\n");
        break; 
    }
    fclose(outputFile);

    printf("Successfully processed %s -> %s\n", inputFilename, outputFilename);
}

// Centers the signal around zero.
void removeDCOffset(double* signal, int numSamples) {
    double sum = 0.0; for (int i = 0; i < numSamples; i++) sum += signal[i];
    double mean = sum / numSamples; for (int i = 0; i < numSamples; i++) signal[i] -= mean;
}

// Scales signal amplitude to [-1, 1].
void normalizeSignal(double* signal, int numSamples) {
    double maxAbsValue = 0.0; for (int i = 0; i < numSamples; i++) if (fabs(signal[i]) > maxAbsValue) maxAbsValue = fabs(signal[i]);
    if (maxAbsValue > 0) for (int i = 0; i < numSamples; i++) signal[i] /= maxAbsValue;
}

// Locates the stable (loudest) part of the vowel.
int findSteadyFrames(const double* signal, int numSamples, int* startFrameIndex) {
    int numFrames = numSamples / FRAME_SIZE; if (numFrames < NUM_STEADY_FRAMES) return 0;
    double maxEnergy = -1.0; int bestStartIndex = -1;
    for (int i = 0; i <= numFrames - NUM_STEADY_FRAMES; i++) {
        double currentBlockEnergy = 0.0;
        for (int j = 0; j < NUM_STEADY_FRAMES; j++) {
            int frameStart = (i + j) * FRAME_SIZE;
            for (int k = 0; k < FRAME_SIZE; k++) currentBlockEnergy += signal[frameStart + k] * signal[frameStart + k];
        }
        if (currentBlockEnergy > maxEnergy) {
            maxEnergy = currentBlockEnergy; bestStartIndex = i * FRAME_SIZE;
        }
    }
    if (bestStartIndex != -1) {*startFrameIndex = bestStartIndex; return 1;} return 0;
}

// Calculates autocorrelation (R values).
void calculateAutocorrelation(const double* frame, double* R) {
    for (int k = 0; k <= LPC_ORDER; k++) {
        R[k] = 0.0; for (int n = 0; n < FRAME_SIZE - k; n++) R[k] += frame[n] * frame[n + k];
    }
}

// Calculates LPC coefficients (A values) from R.
void levinsonDurbin(const double* R, double* A, double* E) {
    double k[LPC_ORDER + 1]; double alpha[LPC_ORDER + 1][LPC_ORDER + 1] = {0};
    if (R[0] == 0.0) {for(int i = 0; i <= LPC_ORDER; ++i) A[i] = 0.0; *E = 0.0; return;}
    *E = R[0];
    for (int i = 1; i <= LPC_ORDER; i++) {
        double sum = 0.0; for (int j = 1; j < i; j++) sum += alpha[i - 1][j] * R[i - j];
        k[i] = (R[i] - sum) / *E; alpha[i][i] = k[i];
        for (int j = 1; j < i; j++) alpha[i][j] = alpha[i - 1][j] - k[i] * alpha[i - 1][i - j];
        *E = (1 - k[i] * k[i]) * *E;
    }
    for (int i = 1; i <= LPC_ORDER; i++) A[i] = alpha[LPC_ORDER][i];
}