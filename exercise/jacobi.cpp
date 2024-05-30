#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

// Function to perform parallel Jacobi iteration
void parallelJacobi(int n, const std::vector<double>& A, const std::vector<double>& b, std::vector<double>& x, int max_iters, double tol) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    assert(n % size == 0); // Ensure n is divisible by size

    int local_n = n / size; // Number of rows per process
    std::vector<double> local_A(local_n * n); // Local part of A
    std::vector<double> local_b(local_n); // Local part of b
    std::vector<double> local_x(local_n); // Local part of x

    // Scatter A and b to local processes
    MPI_Scatter(A.data(), local_n * n, MPI_DOUBLE, local_A.data(), local_n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b.data(), local_n, MPI_DOUBLE, local_b.data(), local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Initialize local x to 0
    std::fill(local_x.begin(), local_x.end(), 0.0);

    std::vector<double> x_old(n, 0.0);

    for (int k = 0; k < max_iters; ++k) {
        // Gather x from all processes
        MPI_Allgather(local_x.data(), local_n, MPI_DOUBLE, x.data(), local_n, MPI_DOUBLE, MPI_COMM_WORLD);

        // Jacobi iteration for local part
        for (int i = 0; i < local_n; ++i) {
            double sigma = 0.0;
            int global_i = rank * local_n + i;
            for (int j = 0; j < n; ++j) {
                if (j != global_i) {
                    sigma += local_A[i * n + j] * x[j];
                }
            }
            local_x[i] = (local_b[i] - sigma) / local_A[i * n + global_i];
        }

        // Check convergence
        double local_diff = 0.0;
        for (int i = 0; i < local_n; ++i) {
            local_diff += (local_x[i] - x[rank * local_n + i]) * (local_x[i] - x[rank * local_n + i]);
        }
        double global_diff = 0.0;
        MPI_Allreduce(&local_diff, &global_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (sqrt(global_diff) < tol) {
            if (rank == 0) {
                std::cout << "Converged after " << k + 1 << " iterations." << std::endl;
            }
            break;
        }

        // Update x_old
        x_old = x;
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int n = 8; // Size of the system
    const int max_iters = 1000; // Maximum number of iterations
    const double tol = 1e-6; // Convergence tolerance

    // Initialize matrix A and vector b
    std::vector<double> A(n * n, 0);
    std::vector<double> b(n, 0);
    std::vector<double> x(n, 0);

    if (rank == 0) {
        // Example system (diagonally dominant for convergence)
        for (int i = 0; i < n; ++i) {
            b[i] = i + 1;
            x[i] = 0;
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    A[i * n + j] = 4;
                } else if (std::abs(i - j) == 1) {
                    A[i * n + j] = 1;
                } else {
                    A[i * n + j] = 0;
                }
            }
        }
    }

    // Perform parallel Jacobi iteration
    parallelJacobi(n, A, b, x, max_iters, tol);

    // Print the result
    if (rank == 0) {
        std::cout << "Result:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    MPI_Finalize();
    return 0;
}
