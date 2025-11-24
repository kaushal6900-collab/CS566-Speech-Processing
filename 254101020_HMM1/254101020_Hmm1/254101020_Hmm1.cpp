#include "StdAfx.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>    // For setting precision
#include <cfloat>     // For DBL_MAX
#include "StdAfx.h" // For Visual Studio 2010

// --- Global HMM Parameters ---
const int N = 5;  // Number of states
const int M = 32; // Number of distinct observation symbols
const int MAX_ITERATIONS = 100; // Max iterations for Baum-Welch
const long double PROB_FLOOR = 1e-30; // To prevent zero probabilities

// Use long double for better precision
using namespace std;
typedef long double ldouble;

// --- Data Structures ---
vector<vector<ldouble>> A(N, vector<ldouble>(N));
vector<vector<ldouble>> B(N, vector<ldouble>(M));
vector<ldouble> pi(N);
vector<int> O;
int T; // Length of observation sequence

vector<vector<ldouble>> alpha;
vector<vector<ldouble>> beta;
vector<vector<ldouble>> gamma;
vector<vector<vector<ldouble>>> xi;
vector<ldouble> c; // Scaling factors

/**
 * Loads the HMM model (A, B, pi) from the "Initial_Model.txt" file.
 */
void loadModel(const string& filename) {
    ifstream file(filename);
    string line;
    if (!file.is_open()) {
        cerr << "Error: Could not open " << filename << endl;
        exit(1);
    }
    while (getline(file, line) && line.find("A Matrix") == string::npos);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) file >> A[i][j];
    }
    while (getline(file, line) && line.find("B Matrix") == string::npos);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) file >> B[i][j];
    }
    while (getline(file, line) && line.find("PI Matrix") == string::npos);
    for (int i = 0; i < N; ++i) file >> pi[i];
    file.close();
}

/**
 * Loads the observation sequence from "HMM_OBSERVATION_SEQUENCE_1.txt".
 */
void loadObservations(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open " << filename << endl;
        exit(1);
    }
    int obs;
    while (file >> obs) {
        O.push_back(obs - 1); // Convert to 0-indexed
    }
    T = O.size();
    file.close();

    // Resize dynamic structures based on T
    alpha.assign(T, vector<ldouble>(N));
    beta.assign(T, vector<ldouble>(N));
    gamma.assign(T, vector<ldouble>(N));
    xi.assign(T, vector<vector<ldouble>>(N, vector<ldouble>(N)));
    c.assign(T, 0.0);
}

/**
 * Writes the final A matrix to "HMM_AIJ_FINAL.txt"
 */
void writeAMatrix(const string& filename) {
    ofstream outFile(filename);
    outFile << fixed << setprecision(10) << scientific;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            outFile << A[i][j] << (j == N - 1 ? "" : "\t");
        }
        outFile << endl;
    }
    outFile.close();
}

/**
 * Writes the final B matrix to "HMM_BJK_FINAL.txt"
 */
void writeBMatrix(const string& filename) {
    ofstream outFile(filename);
    outFile << fixed << setprecision(10) << scientific;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            outFile << B[i][j] << (j == M - 1 ? "" : "\t");
        }
        outFile << endl;
    }
    outFile.close();
}

/**
 * Problem 1: Forward Algorithm with Scaling
 * Computes alpha and scaling factors c. Returns log-probability.
 */
ldouble forward() {
    // 1. Initialization (t = 0)
    c[0] = 0.0;
    for (int i = 0; i < N; ++i) {
        alpha[0][i] = pi[i] * B[i][O[0]];
        c[0] += alpha[0][i];
    }
    // Scale alpha[0]
    for (int i = 0; i < N; ++i) {
        alpha[0][i] /= c[0];
    }

    // 2. Recursion (t = 1 to T-1)
    for (int t = 1; t < T; ++t) {
        c[t] = 0.0;
        for (int j = 0; j < N; ++j) {
            ldouble sum = 0.0;
            for (int i = 0; i < N; ++i) {
                sum += alpha[t - 1][i] * A[i][j];
            }
            alpha[t][j] = sum * B[j][O[t]];
            c[t] += alpha[t][j];
        }
        // Scale alpha[t]
        for (int j = 0; j < N; ++j) {
            alpha[t][j] /= c[t];
        }
    }

    // 3. Termination (Compute Log-Probability)
    ldouble logProb = 0.0;
    for (int t = 0; t < T; ++t) {
        logProb += log(c[t]);
    }
    return logProb;
}

/**
 * Backward Algorithm with Scaling
 * Computes beta, using scaling factors c from forward pass.
 */
void backward() {
    // 1. Initialization (t = T-1)
    for (int i = 0; i < N; ++i) {
        beta[T - 1][i] = 1.0; // Already scaled by c[T-1]
    }

    // 2. Recursion (t = T-2 down to 0)
    for (int t = T - 2; t >= 0; --t) {
        for (int i = 0; i < N; ++i) {
            ldouble sum = 0.0;
            for (int j = 0; j < N; ++j) {
                sum += A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
            }
            beta[t][i] = sum / c[t]; // Scale beta
        }
    }
}

/**
 * E-Step: Compute gamma and xi
 */
void compute_gamma_xi() {
    for (int t = 0; t < T - 1; ++t) {
        ldouble denom = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                denom += alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
            }
        }
        
        for (int i = 0; i < N; ++i) {
            gamma[t][i] = 0.0;
            for (int j = 0; j < N; ++j) {
                xi[t][i][j] = (alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j]) / denom;
                gamma[t][i] += xi[t][i][j];
            }
        }
    }
    // Special case for gamma at t = T-1
    ldouble denom = 0.0;
    for (int i = 0; i < N; i++){
        denom += alpha[T-1][i];
    }
    for(int i = 0; i < N; i++){
        gamma[T-1][i] = alpha[T-1][i] / denom;
    }
}

/**
 * M-Step: Re-estimate model parameters pi, A, B
 */
void re_estimate() {
    // 1. Re-estimate pi
    for (int i = 0; i < N; ++i) {
        pi[i] = gamma[0][i];
    }

    // 2. Re-estimate A
    for (int i = 0; i < N; ++i) {
        ldouble sum_gamma_t = 0.0;
        for (int t = 0; t < T - 1; ++t) {
            sum_gamma_t += gamma[t][i];
        }

        for (int j = 0; j < N; ++j) {
            ldouble sum_xi_t = 0.0;
            for (int t = 0; t < T - 1; ++t) {
                sum_xi_t += xi[t][i][j];
            }
            A[i][j] = sum_xi_t / sum_gamma_t;
        }
    }

    // 3. Re-estimate B
    for (int j = 0; j < N; ++j) {
        ldouble sum_gamma_t = 0.0;
        for (int t = 0; t < T; ++t) {
            sum_gamma_t += gamma[t][j];
        }

        for (int k = 0; k < M; ++k) {
            ldouble sum_gamma_t_k = 0.0;
            for (int t = 0; t < T; ++t) {
                if (O[t] == k) {
                    sum_gamma_t_k += gamma[t][j];
                }
            }
            B[j][k] = sum_gamma_t_k / sum_gamma_t;
            
            // Apply flooring to avoid zero probabilities
            if (B[j][k] < PROB_FLOOR) {
                B[j][k] = PROB_FLOOR;
            }
        }
        // Re-normalize B row to sum to 1 after flooring
        ldouble rowSum = 0.0;
        for(int k=0; k < M; ++k) rowSum += B[j][k];
        for(int k=0; k < M; ++k) B[j][k] /= rowSum;
    }
}

/**
 * Main function to run the Baum-Welch Algorithm
 */
int main() {
    // 1. Load initial model and observations
    loadModel("Initial_Model.txt");
    loadObservations("HMM_OBSERVATION_SEQUENCE_1.txt");

    cout << "Starting Baum-Welch Training..." << endl;
    cout << "T = " << T << ", N = " << N << ", M = " << M << endl;

    ldouble oldLogProb = -DBL_MAX;
    ldouble newLogProb = 0.0;

    // 2. Run iterative training
    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        
        // E-Step
        newLogProb = forward();
        backward();
        compute_gamma_xi();
        
        // M-Step
        re_estimate();

        // Check for convergence
        cout << "Iteration " << iter + 1 << ": Log-Probability = " << fixed << setprecision(20) << newLogProb << endl;
        if (newLogProb <= oldLogProb) {
            cout << "Convergence reached." << endl;
            break;
        }
        oldLogProb = newLogProb;
    }

    // 3. Write final trained model to files
    writeAMatrix("HMM_AIJ_FINAL.txt");
    writeBMatrix("HMM_BJK_FINAL.txt");

    cout << "------------------------------------------" << endl;
    cout << "Training complete." << endl;
    cout << "Final model saved to:" << endl;
    cout << "  - HMM_AIJ_FINAL.txt" << endl;
    cout << "  - HMM_BJK_FINAL.txt" << endl;
    
    // For Visual Studio, to prevent the console from closing
    cout << "Press any key to continue . . ." << endl;
    cin.get();

    return 0;
}