#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdlib>

class DavidsonJacobi {
public:
    DavidsonJacobi(int maxIterations, double tol)
        : maxIterations(maxIterations), tol(tol) {}

    void solve(const std::vector<std::vector<double>>& A, int numEigenvalues,
               std::vector<double>& eigenvalues, std::vector<std::vector<double>>& eigenvectors) {
        int n = A.size();
        eigenvalues.resize(numEigenvalues);
        eigenvectors.resize(n, std::vector<double>(numEigenvalues, 0.0));

        // Initialize V with random values
        std::vector<std::vector<double>> V(n, std::vector<double>(numEigenvalues));
        for (int i = 0; i < numEigenvalues; ++i) {
            for (int j = 0; j < n; ++j) {
                V[j][i] = static_cast<double>(rand()) / RAND_MAX;  // Random initialization
            }
            normalizeColumn(V, i);
        }

        for (int iter = 0; iter < maxIterations; ++iter) {
            // Compute AV
            std::vector<std::vector<double>> AV(n, std::vector<double>(numEigenvalues, 0.0));
            matrixVectorMultiply(A, V, AV);

            // Build the projected matrix H
            std::vector<std::vector<double>> H(numEigenvalues, std::vector<double>(numEigenvalues, 0.0));
            for (int i = 0; i < numEigenvalues; ++i) {
                for (int j = 0; j < numEigenvalues; ++j) {
                    H[i][j] = dotProduct(V, AV, i, j);
                }
            }

            // Eigen decomposition of H (simple method, can be improved)
            std::vector<double> lambda(numEigenvalues);
            std::vector<std::vector<double>> Q(numEigenvalues, std::vector<double>(numEigenvalues));
            eigenDecomposition(H, lambda, Q);

            // Update eigenvalues and eigenvectors
            for (int i = 0; i < numEigenvalues; ++i) {
                eigenvalues[i] = lambda[i];
                for (int j = 0; j < n; ++j) {
                    eigenvectors[j][i] = 0.0;
                    for (int k = 0; k < numEigenvalues; ++k) {
                        eigenvectors[j][i] += V[j][k] * Q[k][i];
                    }
                }
            }

            // Update the basis V
            for (int i = 0; i < numEigenvalues; ++i) {
                for (int j = 0; j < n; ++j) {
                    V[j][i] = 0.0;
                    for (int k = 0; k < numEigenvalues; ++k) {
                        V[j][i] += AV[j][k] * Q[k][i];
                    }
                }
                normalizeColumn(V, i);
            }

            // Check convergence
            if (checkConvergence(V, AV, tol)) {
                break;
            }
        }
    }

private:
    int maxIterations;
    double tol;

    void matrixVectorMultiply(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& V, std::vector<std::vector<double>>& AV) {
        int n = A.size();
        int m = V[0].size();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                AV[i][j] = 0.0;
                for (int k = 0; k < n; ++k) {
                    AV[i][j] += A[i][k] * V[k][j];
                }
            }
        }
    }

    double dotProduct(const std::vector<std::vector<double>>& V, const std::vector<std::vector<double>>& AV, int col1, int col2) {
        double sum = 0.0;
        for (int i = 0; i < V.size(); ++i) {
            sum += V[i][col1] * AV[i][col2];
        }
        return sum;
    }

void normalizeColumn(std::vector<std::vector<double>>& V, int col) {
        double norm = 0.0;
        for (int i = 0; i < V.size(); ++i) {
            norm += V[i][col] * V[i][col];
        }
        norm = std::sqrt(norm);
        for (int i = 0; i < V.size(); ++i) {
            V[i][col] /= norm;
        }
    }

    void eigenDecomposition(const std::vector<std::vector<double>>& H, std::vector<double>& eigenvalues, std::vector<std::vector<double>>& eigenvectors) {
        // This is a placeholder for the eigen decomposition of a small symmetric matrix
        // For demonstration, we will just assign values for simplicity
        for (size_t i = 0; i < eigenvalues.size(); ++i) {
            eigenvalues[i] = H[i][i]; // Diagonal as eigenvalues (not correct in general)
            for (size_t j = 0; j < eigenvectors.size(); ++j) {
                eigenvectors[j][i] = (i == j) ? 1.0 : 0.0; // Identity as eigenvectors (not correct in general)
            }
        }
    }

    bool checkConvergence(const std::vector<std::vector<double>>& V, const std::vector<std::vector<double>>& AV, double tol) {
        double norm = 0.0;
        for (size_t i = 0; i < V.size(); ++i) {
            for (size_t j = 0; j < V[0].size(); ++j) {
                double diff = AV[i][j] - dotProduct(V, V, j, j);
                norm += diff * diff;
            }
        }
        return std::sqrt(norm) < tol;
    }
};
