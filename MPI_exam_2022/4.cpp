/* 在每个进程中，都顺序的向所有进程发送本进程的实型数据，之后顺序的接受来自所有进程
的实型数据，构成一个 n × n 的阵列 (n 为进程数，指每个进程中都有一个长度为 n 的
实型数组，存储来自这 n 个进程的数据) */

#include <iostream>
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int block_size = 2;  // 假设每个块包含2个double
    double send_data[block_size];
    double recv_data[block_size * size];

    // 初始化发送数据
    for (int i = 0; i < block_size; ++i) {
        send_data[i] = rank * 10 + i;
    }

    // 使用 MPI_Send 和 MPI_Recv 实现 MPI_Allgather
    for (int i = 0; i < size; ++i) {
        if (rank == i) {
            // 当前进程发送数据给所有其他进程
            for (int j = 0; j < size; ++j) {
                if (j != i) {
                    MPI_Send(send_data, block_size, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
                }
            }
        } else {
            // 从发送进程接收数据
            MPI_Recv(&recv_data[i * block_size], block_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // 当前进程自身的数据复制到接收缓冲区
    for (int i = 0; i < block_size; ++i) {
        recv_data[rank * block_size + i] = send_data[i];
    }

    // 打印所有接收到的数据
    std::cout << "Process " << rank << " received data: ";
    for (int i = 0; i < block_size * size; ++i) {
        std::cout << recv_data[i] << " ";
    }
    std::cout << std::endl;

    MPI_Finalize();
    return 0;
}
