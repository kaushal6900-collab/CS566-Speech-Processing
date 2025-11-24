#ifndef HMM_H
#define HMM_H

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

const int NUM_STATES = 5;
// CODEBOOK_SIZE is defined in Codebook.h, but we need it here. 
// To avoid circular dependency or redefinition, we can redefine or include Codebook.h.
// Or just use 64.
const int HMM_CODEBOOK_SIZE = 64;

class HMM {
public:
    long double A[NUM_STATES][NUM_STATES];
    long double B[NUM_STATES][HMM_CODEBOOK_SIZE];
    long double Pi[NUM_STATES];

    HMM() {
        Initialize();
    }

    void Initialize() {
        // Initialize Pi: [1, 0, 0, 0, 0]
        for (int i = 0; i < NUM_STATES; ++i) {
            Pi[i] = (i == 0) ? 1.0 : 0.0;
        }

        // Initialize A: Bakis model (left-to-right)
        // Simple initialization: 0.5 to self, 0.5 to next
        for (int i = 0; i < NUM_STATES; ++i) {
            for (int j = 0; j < NUM_STATES; ++j) {
                A[i][j] = 0.0;
            }
        }

        for (int i = 0; i < NUM_STATES; ++i) {
            if (i < NUM_STATES - 1) {
                A[i][i] = 0.5;
                A[i][i + 1] = 0.5;
            } else {
                A[i][i] = 1.0; // Last state absorbs
            }
        }

        // Initialize B: Uniform
        long double uniform_prob = 1.0 / HMM_CODEBOOK_SIZE;
        for (int i = 0; i < NUM_STATES; ++i) {
            for (int k = 0; k < HMM_CODEBOOK_SIZE; ++k) {
                B[i][k] = uniform_prob;
            }
        }
    }

    // Returns log-likelihood
    long double Forward(const vector<int>& O, vector<vector<long double>>& alpha, vector<long double>& scale) {
        int T = O.size();
        alpha.assign(T, vector<long double>(NUM_STATES, 0.0));
        scale.assign(T, 0.0);

        // Initialization
        for (int i = 0; i < NUM_STATES; ++i) {
            alpha[0][i] = Pi[i] * B[i][O[0]];
            scale[0] += alpha[0][i];
        }

        // Scaling t=0
        if (scale[0] <= 0) scale[0] = 1e-30;
        for (int i = 0; i < NUM_STATES; ++i) {
            alpha[0][i] /= scale[0];
        }

        // Induction
        for (int t = 1; t < T; ++t) {
            scale[t] = 0.0;
            for (int j = 0; j < NUM_STATES; ++j) {
                long double sum = 0.0;
                for (int i = 0; i < NUM_STATES; ++i) {
                    sum += alpha[t - 1][i] * A[i][j];
                }
                alpha[t][j] = sum * B[j][O[t]];
                scale[t] += alpha[t][j];
            }
            // Scaling
            if (scale[t] <= 0) scale[t] = 1e-30;
            for (int j = 0; j < NUM_STATES; ++j) {
                alpha[t][j] /= scale[t];
            }
        }

        // Termination (Log Likelihood)
        long double logLikelihood = 0.0;
        for (int t = 0; t < T; ++t) {
            logLikelihood += log(scale[t]);
        }
        return logLikelihood;
    }

    void Backward(const vector<int>& O, vector<vector<long double>>& beta, const vector<long double>& scale) {
        int T = O.size();
        beta.assign(T, vector<long double>(NUM_STATES, 0.0));

        // Initialization
        for (int i = 0; i < NUM_STATES; ++i) {
            beta[T - 1][i] = 1.0 / scale[T - 1]; 
        }

        // Induction
        for (int t = T - 2; t >= 0; --t) {
            for (int i = 0; i < NUM_STATES; ++i) {
                long double sum = 0.0;
                for (int j = 0; j < NUM_STATES; ++j) {
                    sum += A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
                }
                beta[t][i] = sum / scale[t]; 
            }
        }
    }

    void BaumWelch(const vector<vector<int>>& observations) {
        int num_iter = 30;
        int num_obs = observations.size();

        for (int iter = 0; iter < num_iter; ++iter) {
            long double expected_A_num[NUM_STATES][NUM_STATES] = {0};
            long double expected_A_den[NUM_STATES] = {0};
            long double expected_B_num[NUM_STATES][HMM_CODEBOOK_SIZE] = {0};
            long double expected_B_den[NUM_STATES] = {0};
            long double pi_sum[NUM_STATES] = {0};

            for (int r = 0; r < num_obs; ++r) {
                const vector<int>& O = observations[r];
                int T = O.size();
                vector<vector<long double>> alpha, beta;
                vector<long double> scale;

                Forward(O, alpha, scale);
                Backward(O, beta, scale);

                vector<vector<long double>> gamma(T, vector<long double>(NUM_STATES));
                for(int t=0; t<T; ++t) {
                    long double denom = 0.0;
                    for(int i=0; i<NUM_STATES; ++i) {
                        denom += alpha[t][i] * beta[t][i];
                    }
                    for(int i=0; i<NUM_STATES; ++i) {
                        gamma[t][i] = (alpha[t][i] * beta[t][i]) / denom;
                    }
                }

                // Accumulate Pi
                for(int i=0; i<NUM_STATES; ++i) {
                    pi_sum[i] += gamma[0][i];
                }

                // Accumulate A
                for(int t=0; t<T-1; ++t) {
                    long double denom = 0.0;
                    for(int i=0; i<NUM_STATES; ++i) {
                        for(int j=0; j<NUM_STATES; ++j) {
                            denom += alpha[t][i] * A[i][j] * B[j][O[t+1]] * beta[t+1][j];
                        }
                    }
                    
                    for(int i=0; i<NUM_STATES; ++i) {
                        expected_A_den[i] += gamma[t][i];
                        for(int j=0; j<NUM_STATES; ++j) {
                            long double xi = (alpha[t][i] * A[i][j] * B[j][O[t+1]] * beta[t+1][j]) / denom;
                            expected_A_num[i][j] += xi;
                        }
                    }
                }

                // Accumulate B
                for(int t=0; t<T; ++t) {
                    for(int i=0; i<NUM_STATES; ++i) {
                        expected_B_den[i] += gamma[t][i];
                        expected_B_num[i][O[t]] += gamma[t][i];
                    }
                }
            }

            // Update Parameters
            for(int i=0; i<NUM_STATES; ++i) {
                Pi[i] = pi_sum[i] / num_obs;
            }

            for(int i=0; i<NUM_STATES; ++i) {
                for(int j=0; j<NUM_STATES; ++j) {
                    if (expected_A_den[i] > 0)
                        A[i][j] = expected_A_num[i][j] / expected_A_den[i];
                    else
                        A[i][j] = (i==j) ? 1.0 : 0.0;
                }
            }

            for(int i=0; i<NUM_STATES; ++i) {
                for(int k=0; k<HMM_CODEBOOK_SIZE; ++k) {
                    if (expected_B_den[i] > 0) {
                        B[i][k] = expected_B_num[i][k] / expected_B_den[i];
                        if (B[i][k] < 1e-30) B[i][k] = 1e-30;
                    } else {
                        B[i][k] = 1.0 / HMM_CODEBOOK_SIZE;
                    }
                }
            }
        }
    }

    long double GetProbability(const vector<int>& O) {
        vector<vector<long double>> alpha;
        vector<long double> scale;
        return Forward(O, alpha, scale);
    }

    void Save(string filename) {
        ofstream out(filename);
        if (!out) {
            cerr << "Error saving model to " << filename << endl;
            return;
        }
        out << scientific << setprecision(10);
        for (int i = 0; i < NUM_STATES; ++i) out << Pi[i] << " ";
        out << endl;
        for (int i = 0; i < NUM_STATES; ++i) {
            for (int j = 0; j < NUM_STATES; ++j) out << A[i][j] << " ";
            out << endl;
        }
        for (int i = 0; i < NUM_STATES; ++i) {
            for (int k = 0; k < HMM_CODEBOOK_SIZE; ++k) out << B[i][k] << " ";
            out << endl;
        }
        out.close();
    }

    void Load(string filename) {
        ifstream in(filename);
        if (!in) {
            cerr << "Error loading model from " << filename << endl;
            return;
        }
        for (int i = 0; i < NUM_STATES; ++i) in >> Pi[i];
        for (int i = 0; i < NUM_STATES; ++i) {
            for (int j = 0; j < NUM_STATES; ++j) in >> A[i][j];
        }
        for (int i = 0; i < NUM_STATES; ++i) {
            for (int k = 0; k < HMM_CODEBOOK_SIZE; ++k) in >> B[i][k];
        }
        in.close();
    }
};

#endif