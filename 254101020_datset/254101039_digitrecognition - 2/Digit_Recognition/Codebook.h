#pragma once
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

// Constants
#define CODEBOOK_SIZE 64
#define EPSILON 0.03
#define DELTA 0.00001

class Codebook {
public:
    // Stores the 64 centroids (each dim 12)
    vector<vector<double> > centroids;

    // Euclidean distance (squared)
    double ComputeDistance(const vector<double>& v1, const vector<double>& v2) {
        double dist = 0.0;
        for (size_t i = 0; i < v1.size(); ++i) {
            double diff = v1[i] - v2[i];
            dist += diff * diff;
        }
        return dist;
    }

    // Implements LBG Algorithm
    void Train(const vector<vector<double> >& universe) {
        if (universe.empty()) return;
        size_t dimension = universe[0].size();

        // 1. Start with 1 centroid (mean of universe)
        vector<double> global_mean(dimension, 0.0);
        for (size_t k = 0; k < universe.size(); ++k) {
            const vector<double>& vec = universe[k];
            for (size_t i = 0; i < dimension; ++i) {
                global_mean[i] += vec[i];
            }
        }
        for (size_t i = 0; i < dimension; ++i) {
            global_mean[i] /= universe.size();
        }

        centroids.clear();
        centroids.push_back(global_mean);

        // 2. Loop: Split centroids -> 2, 4, 8... 64
        while (centroids.size() < CODEBOOK_SIZE) {
            // Split existing centroids
            vector<vector<double> > next_centroids;
            next_centroids.reserve(centroids.size() * 2);
            
            for (size_t k = 0; k < centroids.size(); ++k) {
                vector<double> c1 = centroids[k];
                vector<double> c2 = centroids[k];
                for (size_t i = 0; i < dimension; ++i) {
                    c1[i] *= (1.0 + EPSILON);
                    c2[i] *= (1.0 - EPSILON);
                }
                next_centroids.push_back(c1);
                next_centroids.push_back(c2);
            }
            centroids = next_centroids;

            // Inner Loop (K-means)
            double prev_distortion = numeric_limits<double>::max();
            
            while (true) {
                vector<vector<double> > sums(centroids.size(), vector<double>(dimension, 0.0));
                vector<int> counts(centroids.size(), 0);
                double current_distortion = 0.0;

                // Assign vectors to nearest centroid
                for (size_t k = 0; k < universe.size(); ++k) {
                    const vector<double>& vec = universe[k];
                    int nearest_idx = -1;
                    double min_dist = numeric_limits<double>::max();

                    for (size_t i = 0; i < centroids.size(); ++i) {
                        double d = ComputeDistance(vec, centroids[i]);
                        if (d < min_dist) {
                            min_dist = d;
                            nearest_idx = i;
                        }
                    }

                    current_distortion += min_dist;
                    counts[nearest_idx]++;
                    for (size_t i = 0; i < dimension; ++i) {
                        sums[nearest_idx][i] += vec[i];
                    }
                }

                // Calculate average distortion
                current_distortion /= universe.size();

                // Check distortion change
                if (prev_distortion != numeric_limits<double>::max()) {
                    double improvement = (prev_distortion - current_distortion) / prev_distortion;
                    if (improvement < DELTA) {
                        break; // Converged
                    }
                }
                prev_distortion = current_distortion;

                // Update centroids
                for (size_t i = 0; i < centroids.size(); ++i) {
                    if (counts[i] > 0) {
                        for (size_t j = 0; j < dimension; ++j) {
                            centroids[i][j] = sums[i][j] / counts[i];
                        }
                    }
                }
            }
        }
    }

    // Returns index (0-63) of nearest centroid
    int Quantize(const vector<double>& feature) {
        int nearest_idx = -1;
        double min_dist = numeric_limits<double>::max();

        for (size_t i = 0; i < centroids.size(); ++i) {
            double d = ComputeDistance(feature, centroids[i]);
            if (d < min_dist) {
                min_dist = d;
                nearest_idx = i;
            }
        }
        return nearest_idx;
    }

    // Save centroids to file
    void Save(string filename) {
        ofstream outfile(filename.c_str());
        if (!outfile.is_open()) {
            cerr << "Error opening file for saving: " << filename << endl;
            return;
        }

        if (!centroids.empty()) {
            outfile << centroids.size() << " " << centroids[0].size() << endl;
            for (size_t k = 0; k < centroids.size(); ++k) {
                const vector<double>& c = centroids[k];
                for (size_t i = 0; i < c.size(); ++i) {
                    outfile << c[i] << (i == c.size() - 1 ? "" : " ");
                }
                outfile << endl;
            }
        }
        outfile.close();
    }

    // Load centroids from file
    void Load(string filename) {
        ifstream infile(filename.c_str());
        if (!infile.is_open()) {
            cerr << "Error opening file for loading: " << filename << endl;
            return;
        }

        int size, dim;
        if (infile >> size >> dim) {
            centroids.resize(size, vector<double>(dim));
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < dim; ++j) {
                    infile >> centroids[i][j];
                }
            }
        }
        infile.close();
    }
};