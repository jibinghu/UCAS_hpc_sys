// 实现矩阵A和向量b按行存储在p个处理器中的Ax+b算法
#include <mpi.h>
#include <vector>
#include <iostream>

void matrix_vector_multiply(int rows, int cols, const std::vector<int>& local_A, const std::vector<int>& x, std::vector<int>& local_y) {
    for (int i = 0; i < rows; ++i) {
        local_y[i] = 0; // Initialize local_y[i] to zero before adding
        for (int j = 0; j < cols; ++j) {
            local_y[i] += local_A[i * cols + j] * x[j];
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int m = 4; // Number of rows per process (assuming each process has the same number of rows)
    int n = m * size; // Total number of rows/columns (assuming a square matrix for simplicity)

    // Each process already has its submatrix of size m x n and subvector of size m
    std::vector<int> local_A(m * n); // Submatrix
    std::vector<int> local_b(m);     // Subvector b
    std::vector<int> local_y(m);     // Resultant subvector

    // Initialize local_A and local_b with some values (for demonstration purposes)
    for (int i = 0; i < m * n; ++i) {
        local_A[i] = rank + 1; // Example values
    }
    for (int i = 0; i < m; ++i) {
        local_b[i] = rank + 1; // Example values
    }

    // Global vector x (only initialized in rank 0 and broadcasted to all processes)
    std::vector<int> x(n);
    if (rank == 0) {
        // Initialize global vector x with some values (for demonstration purposes)
        for (int i = 0; i < n; ++i) {
            x[i] = i + 1; // Example values
        }
    }

    // Broadcast the vector x to all processes
    MPI_Bcast(x.data(), n, MPI_INT, 0, MPI_COMM_WORLD);

    // Perform the local matrix-vector multiplication
    matrix_vector_multiply(m, n, local_A, x, local_y);

    // Add local vector b to the result
    for (int i = 0; i < m; ++i) {
        local_y[i] += local_b[i];
    }

    // Gather all partial results y to get the final result
    std::vector<int> global_y(n);
    MPI_Gather(local_y.data(), m, MPI_INT, global_y.data(), m, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Print the result
        std::cout << "Resulting vector x+ = A*x + b:\n";
        for (int i = 0; i < n; ++i) {
            std::cout << global_y[i] << " ";
        }
        std::cout << std::endl;
    }

    MPI_Finalize();
    return 0;
}
