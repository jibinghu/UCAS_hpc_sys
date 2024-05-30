#include <mpi.h>
#include <iostream>
#include <vector>

void initializeMatrices(int n, int p, std::vector<int>& local_A, std::vector<int>& local_B, std::vector<int>& local_C) {
    int blockSize = n / p;
    for (int i = 0; i < blockSize; ++i) {
        for (int j = 0; j < blockSize; ++j) {
            local_A[i * blockSize + j] = rand() % 10; // Example initialization
            local_B[i * blockSize + j] = rand() % 10;
            local_C[i * blockSize + j] = 0;
        }
    }
}

void matrixMultiply(int n, int p, std::vector<int>& local_A, std::vector<int>& local_B, std::vector<int>& local_C) {
    int blockSize = n / p;
    std::vector<int> temp_A(blockSize * blockSize);
    std::vector<int> temp_B(blockSize * blockSize);

    for (int iter = 0; iter < p; ++iter) {
        // Perform local matrix multiplication and accumulate the result
        for (int i = 0; i < blockSize; ++i) {
            for (int j = 0; j < blockSize; ++j) {
                for (int k = 0; k < blockSize; ++k) {
                    local_C[i * blockSize + j] += local_A[i * blockSize + k] * local_B[k * blockSize + j];
                }
            }
        }

        // Shift local_A left by one position
        MPI_Sendrecv_replace(local_A.data(), blockSize * blockSize, MPI_INT, 
                             (MPI::COMM_WORLD.Get_rank() - 1 + p) % p, 0, 
                             (MPI::COMM_WORLD.Get_rank() + 1) % p, 0, 
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Shift local_B up by one position
        MPI_Sendrecv_replace(local_B.data(), blockSize * blockSize, MPI_INT, 
                             (MPI::COMM_WORLD.Get_rank() - p + p) % p, 0, 
                             (MPI::COMM_WORLD.Get_rank() + p) % p, 0, 
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int n = 4; // Assume n is divisible by p
    const int p = 2; // Number of processes per dimension
    const int blockSize = n / p;

    std::vector<int> local_A(blockSize * blockSize);
    std::vector<int> local_B(blockSize * blockSize);
    std::vector<int> local_C(blockSize * blockSize, 0);

    initializeMatrices(n, p, local_A, local_B, local_C);

    // Cannon's algorithm
    matrixMultiply(n, p, local_A, local_B, local_C);

    // Print result
    if (rank == 0) {
        std::cout << "Result matrix C:" << std::endl;
        for (int i = 0; i < blockSize; ++i) {
            for (int j = 0; j < blockSize; ++j) {
                std::cout << local_C[i * blockSize + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}
