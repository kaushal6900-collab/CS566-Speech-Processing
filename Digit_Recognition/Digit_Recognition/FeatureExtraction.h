#ifndef FEATURE_EXTRACTION_H
#define FEATURE_EXTRACTION_H

#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

// Constants
const int FRAME_SIZE = 320;
const int FRAME_SHIFT = 80;
const int P = 12; // LPC Order
const int Q = 12; // Cepstral Order
const double PI = 3.14159265358979323846;

// Helper function: Apply Hamming Window
inline void ApplyHammingWindow(vector<double>& frame) {
    int N = frame.size();
    for (int i = 0; i < N; ++i) {
        frame[i] *= (0.54 - 0.46 * cos(2 * PI * i / (N - 1)));
    }
}

// Helper function: Compute Autocorrelation
inline void ComputeAutoCorrelation(const vector<double>& frame, vector<double>& R) {
    R.assign(P + 1, 0.0);
    int N = frame.size();
    for (int k = 0; k <= P; ++k) {
        for (int m = 0; m < N - k; ++m) {
            R[k] += frame[m] * frame[m + k];
        }
    }
}

// Helper function: Durbin's Algorithm
// Returns the error energy (gain squared)
inline double DurbinAlgo(const vector<double>& R, vector<double>& alpha) {
    double E = R[0];
    alpha.assign(P + 1, 0.0);
    
    // Temporary vectors for recursion
    vector<double> alpha_prev(P + 1, 0.0);
    
    // Use long double for internal calculations to avoid instability
    long double k; 

    for (int i = 1; i <= P; ++i) {
        long double sum = 0.0;
        for (int j = 1; j <= i - 1; ++j) {
            sum += alpha_prev[j] * R[i - j];
        }
        
        k = (R[i] - sum) / E;
        
        alpha[i] = (double)k;
        
        for (int j = 1; j <= i - 1; ++j) {
            alpha[j] = alpha_prev[j] - (double)k * alpha_prev[i - j];
        }
        
        E = (1 - k * k) * E;
        
        // Update alpha_prev for next iteration
        for(int j = 1; j <= i; ++j) {
            alpha_prev[j] = alpha[j];
        }
    }
    
    return E;
}

// Helper function: Compute Cepstral Coefficients
inline void ComputeCepstral(const vector<double>& alpha, vector<double>& C, double gain) {
    C.assign(Q + 1, 0.0);
    
    // C[0] is typically log energy
    if (gain > 0)
        C[0] = log(gain);
    else
        C[0] = 0;

    for (int m = 1; m <= Q; ++m) {
        double sum = 0.0;
        for (int k = 1; k < m; ++k) {
            sum += (double)k / m * C[k] * alpha[m - k];
        }
        C[m] = alpha[m] + sum;
    }
}

// Helper function: Apply Liftering
inline void ApplyLiftering(vector<double>& C) {
    for (int m = 1; m <= Q; ++m) {
        C[m] *= (1 + (Q / 2.0) * sin(PI * m / Q));
    }
}

// Main function: Compute Features from Samples
inline void ComputeFeatures(const vector<double>& input_samples, vector<vector<double> >& all_features) {
    if (input_samples.empty()) return;

    vector<double> samples = input_samples; // Copy to modify

    // 1. DC Shift
    double mean = 0.0;
    for (size_t i = 0; i < samples.size(); ++i) mean += samples[i];
    mean /= samples.size();
    for (size_t i = 0; i < samples.size(); ++i) samples[i] -= mean;

    // 2. Normalization
    double max_amp = 0.0;
    for (size_t i = 0; i < samples.size(); ++i) {
        if (abs(samples[i]) > max_amp) max_amp = abs(samples[i]);
    }
    if (max_amp > 0) {
        double scale = 5000.0 / max_amp;
        for (size_t i = 0; i < samples.size(); ++i) samples[i] *= scale;
    }

    // 3. Framing and Processing
    int num_samples = samples.size();
    
    for (int start = 0; start + FRAME_SIZE <= num_samples; start += FRAME_SHIFT) {
        vector<double> frame(FRAME_SIZE);
        for (int i = 0; i < FRAME_SIZE; ++i) {
            frame[i] = samples[start + i];
        }

        // Apply Hamming Window
        ApplyHammingWindow(frame);

        // Compute Autocorrelation
        vector<double> R;
        ComputeAutoCorrelation(frame, R);

        // Durbin's Algorithm
        vector<double> alpha;
        double E = DurbinAlgo(R, alpha);

        // Compute Cepstral Coefficients
        vector<double> C;
        ComputeCepstral(alpha, C, E);

        // Apply Liftering
        ApplyLiftering(C);

        // Store features (Storing C[1] to C[Q])
        vector<double> features;
        for(int i = 1; i <= Q; ++i) {
            features.push_back(C[i]);
        }
        all_features.push_back(features);
    }
}

// Helper: Get Feature Vectors from File
inline void GetFeatureVectors(string filePath, vector<vector<double> >& all_features) {
    ifstream infile(filePath.c_str());
    if (!infile.is_open()) {
        cerr << "Error opening file: " << filePath << endl;
        return;
    }

    vector<double> samples;
    double val;
    while (infile >> val) {
        samples.push_back(val);
    }
    infile.close();

    ComputeFeatures(samples, all_features);
}

#endif // FEATURE_EXTRACTION_H