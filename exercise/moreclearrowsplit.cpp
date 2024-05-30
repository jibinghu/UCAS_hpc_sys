// 假设有p pp个进程，A AA是n ∗ n n*nn∗n阶矩阵，n = m ∗ p n=m*pn=m∗p，设计并行计算A x + b Ax+bAx+b的算法，其中A AA是按列分块存放在每个进程中，亦即每个进程中的A AA是n ∗ m n*mn∗m阶矩阵，每个进程中x xx是m mm维向量，b bb存放在第0 00个进程中，最终计算结果存放在x xx中。

// （20分） 并行算法简述；
// （45分）写出只使用MPI_Send和MPI_Recv或者MPI_Sendrecv函数实现并行计算A x + b Ax+bAx+b的子程序，计算结果存放在x xx中；
// （35分）当原始矩阵A的元素a i j = i + j ， x = 1 , − 1 , 1 , − 1 , . . . . . . a_{ij}=i+j，x={1, -1, 1, -1, ......}a 
// ij
// ​
//  =i+j，x=1,−1,1,−1,......，b = x b=xb=x。请写一个主程序对上述子程序进行验证，并给出p = 4 p=4p=4，m = 11 m=11m=11时的计算结果。

#include <mpi.h>
#include <iostream>
#include <vector>

void parallelAxPlusB(int n, int m, const std::vector<int>& local_A, const std::vector<int>& local_x, std::vector<int>& b, std::vector<int>& result) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Step 1: Initialize local result vector
    std::vector<int> local_result(n, 0);

    // Step 2: Compute local A*x
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            local_result[i] += local_A[i * m + j] * local_x[j];
        }
    }

    // Step 3: Send local results to process 0
    if (rank == 0) {
        result = local_result;
        for (int i = 1; i < size; ++i) {
            std::vector<int> temp_result(n);
            MPI_Recv(temp_result.data(), n, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < n; ++j) {
                result[j] += temp_result[j];
            }
        }
        // Add vector b
        for (int i = 0; i < n; ++i) {
            result[i] += b[i];
        }
    } else {
        MPI_Send(local_result.data(), n, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int p = 4;
    const int m = 11;
    const int n = m * p;

    // Initialize matrix A, vector x and vector b
    std::vector<int> local_A(n * m);
    std::vector<int> local_x(m);
    std::vector<int> b(n);
    std::vector<int> result;

    // Fill local_A with values a_ij = i + j
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            local_A[i * m + j] = i + rank * m + j + 1;
        }
    }

    // Fill local_x with alternating 1 and -1
    for (int i = 0; i < m; ++i) {
        local_x[i] = (i % 2 == 0) ? 1 : -1;
    }

    // Fill b with the same values as x in process 0
    if (rank == 0) {
        for (int i = 0; i < n; ++i) {
            b[i] = (i % 2 == 0) ? 1 : -1;
        }
    }

    // Broadcast b to all processes
    MPI_Bcast(b.data(), n, MPI_INT, 0, MPI_COMM_WORLD);

    // Perform parallel computation of A*x + b
    parallelAxPlusB(n, m, local_A, local_x, b, result);

    // Process 0 prints the result
    if (rank == 0) {
        std::cout << "Result:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << result[i] << " ";
        }
        std::cout << std::endl;
    }

    MPI_Finalize();
    return 0;
}
