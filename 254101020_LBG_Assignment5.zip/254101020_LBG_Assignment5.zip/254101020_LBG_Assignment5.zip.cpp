#include "stdafx.h" 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


#define K_FINAL 8
#define DELTA 0.0001
#define EPSILON 0.03
#define VECTOR_DIM 12
#define MAX_M 1100 
#define FILE_NAME "Universe.csv"

// Specified Tokhura Distance Weights
const double TOKHURA_WEIGHTS[VECTOR_DIM] = {
    1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0
};


double universe[MAX_M][VECTOR_DIM];
int M = 0; 



int read_universe(const char *filename) {
    FILE *file;
    
    if (fopen_s(&file, filename, "r") != 0) {
        fprintf(stderr, "Error: Could not open file %s. Ensure the file is in the correct directory.\n", filename);
        return 0;
    }

    M = 0;
    while (M < MAX_M) {
        int count = 0;
        
        for (int i = 0; i < VECTOR_DIM; i++) {
            if (fscanf(file, "%lf%*[, \n]", &universe[M][i]) == 1) {
                count++;
            } else {
                break;
            }
        }
        
        if (count == VECTOR_DIM) {
            M++;
        } else if (count == 0) {
            break; // Reached EOF or an empty line
        } else {
            fprintf(stderr, "Error: Incomplete vector found at line %d.\n", M + 1);
            fclose(file);
            return 0;
        }
    }
    
    fclose(file);
    return 1;
}

double tokhura_distance(const double x[], const double y[]) {
    double distance = 0.0;
    for (int i = 0; i < VECTOR_DIM; i++) {
        distance += TOKHURA_WEIGHTS[i] * (x[i] - y[i]) * (x[i] - y[i]);
    }
    return distance;
}

// Forward declaration of K-Means
double k_means(double codebook[][VECTOR_DIM], int k_current, double (*universe_data)[VECTOR_DIM], int M_size);


// --- K-Means Function (Definition) ---

double k_means(double codebook[][VECTOR_DIM], int k_current, double (*universe_data)[VECTOR_DIM], int M_size) {
    
    // Dynamically allocate memory 
    int *assignments = (int *)malloc(M_size * sizeof(int));
    double (*new_codebook)[VECTOR_DIM] = (double (*)[VECTOR_DIM])malloc(k_current * sizeof(double[VECTOR_DIM]));
    int *counts = (int *)calloc(k_current, sizeof(int));

    if (assignments == NULL || new_codebook == NULL || counts == NULL) {
        fprintf(stderr, "Error: Memory allocation failed in k_means.\n");
        exit(EXIT_FAILURE);
    }
    
    double distortion_prev = DBL_MAX;
    int m = 0;
    
    printf("  Starting K-Means refinement with k=%d...\n", k_current);

    while (1) {
        m++;
        double current_distortion = 0.0;

        // 1. Vector Assignment (Nearest Neighbor Rule)
        for (int i = 0; i < M_size; i++) {
            double min_dist = DBL_MAX;
            int nearest_centroid_index = 0;
            
            for (int j = 0; j < k_current; j++) {
                double dist = tokhura_distance(universe_data[i], codebook[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                    nearest_centroid_index = j;
                }
            }
            
            assignments[i] = nearest_centroid_index;
            current_distortion += min_dist;
        }
        
        // 2. Calculate Distortion

        // 3. Compute and Update Centroids
        for(int j = 0; j < k_current; j++) {
            for(int l = 0; l < VECTOR_DIM; l++) {
                new_codebook[j][l] = 0.0;
            }
            counts[j] = 0;
        }

        for (int i = 0; i < M_size; i++) {
            int j = assignments[i];
            for (int l = 0; l < VECTOR_DIM; l++) {
                new_codebook[j][l] += universe_data[i][l];
            }
            counts[j]++;
        }
        
        for (int j = 0; j < k_current; j++) {
            if (counts[j] > 0) {
                for (int l = 0; l < VECTOR_DIM; l++) {
                    codebook[j][l] = new_codebook[j][l] / counts[j];
                }
            }
            // Empty region: keep the old centroid.
        }

        double avg_distortion = current_distortion / M_size;
        printf("    K-Means Iteration %d: Total Distortion = %.6f, Avg Distortion = %.6f\n", m, current_distortion, avg_distortion);

        // 4. Check Convergence
        if (m > 1 && fabs(distortion_prev - current_distortion) <= DELTA) {
            printf("  K-Means converged after %d iterations (Delta <= %.4f).\n", m, DELTA);
            break;
        }
        
        distortion_prev = current_distortion;
    }

    free(assignments);
    free(new_codebook);
    free(counts);

    return distortion_prev;
}


// --- LBG Algorithm Core ---

double lbg_algorithm(double (*codebook)[VECTOR_DIM]) {
    
    printf("--- Starting LBG Algorithm (Target k=%d) ---\n", K_FINAL);

    // 1. Start with a single vector (centroid of the universe)
    printf("\nInitializing k=1 centroid...\n");
    for (int l = 0; l < VECTOR_DIM; l++) {
        double sum = 0.0;
        for (int i = 0; i < M; i++) {
            sum += universe[i][l];
        }
        codebook[0][l] = sum / M;
    }
    int k_current = 1;
    double final_distortion;
    
    // Initial K-Means refinement
    final_distortion = k_means(codebook, k_current, universe, M);
    
    while (k_current < K_FINAL) {
        // 2. Double the codebook size (Splitting)
        int k_prev = k_current;
        k_current *= 2;
        
        // Split each current centroid
        for (int i = k_prev - 1; i >= 0; i--) {
            // Save the original centroid 
            double Y_orig[VECTOR_DIM];
            for (int l = 0; l < VECTOR_DIM; l++) {
                Y_orig[l] = codebook[i][l];
            }
            
            // New centroid Y+ (overwrites original centroid's location)
            for (int l = 0; l < VECTOR_DIM; l++) {
                codebook[i][l] = Y_orig[l] * (1.0 + EPSILON);
            }

            // New centroid Y- (at location i + k_prev)
            for (int l = 0; l < VECTOR_DIM; l++) {
                codebook[i + k_prev][l] = Y_orig[l] * (1.0 - EPSILON);
            }
        }
        
        printf("\nCodebook doubled from %d to %d. Now performing K-Means refinement.\n", k_prev, k_current);
        
        // 3. Apply K-Means to refine the new centroids
        final_distortion = k_means(codebook, k_current, universe, M);
    }
    
    // Print final codebook
    printf("\nFinal Codebook (k=%d, %d values per vector):\n", K_FINAL, VECTOR_DIM);
    for (int i = 0; i < K_FINAL; i++) {
        printf("Centroid %d: {", i + 1);
        for (int l = 0; l < VECTOR_DIM; l++) {
            printf("%.6f", codebook[i][l]);
            if (l < VECTOR_DIM - 1) {
                printf(", ");
            }
        }
        printf("}\n");
    }

    return final_distortion;
}

// --- Main Execution ---
int main() {
    if (!read_universe(FILE_NAME)) {
        return 1;
    }
    
    printf("Data Loaded: Universe Size (M) = %d, Vector Dimension = %d\n", M, VECTOR_DIM);

    double (*codebook)[VECTOR_DIM] = (double (*)[VECTOR_DIM])malloc(K_FINAL * sizeof(double[VECTOR_DIM]));
    if (codebook == NULL) {
        fprintf(stderr, "Error: Memory allocation failed in main.\n");
        return 1;
    }

    double final_distortion = lbg_algorithm(codebook);
    double final_avg_distortion = final_distortion / M;

    // Output Results
   
    printf("             LBG FINAL RESULTS (k=%d)\n", K_FINAL);
  
    printf("Final Total Distortion (dist): %.6f\n", final_distortion);
    printf("Final Average Distortion (avg_dist): %.6f\n", final_avg_distortion);
   
    
    free(codebook);
    return 0;
}