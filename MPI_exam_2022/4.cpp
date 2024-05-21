#include <iostream>
#include <mpi.h>

// 自定义 MPI_Allgather 实现
void mpi_allgatherSelfImpl(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                           void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
    int rank, nproc;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);

    // 当前进程发送数据给所有其他进程，包括自己
    for (int j = 0; j < nproc; ++j) {
        MPI_Send(sendbuf, sendcount, sendtype, j, 0, comm);
    }

    // 每个进程接收来自所有进程的数据
    for (int i = 0; i < nproc; ++i) {
        MPI_Status status;
        MPI_Recv((char *)recvbuf + i * recvcount * sizeof(float), recvcount, recvtype, i, 0, comm, &status);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int dataLength = 1; // 每个进程发送和接收的数据块大小
    float sendData = rank + 1.0; // 每个进程发送的数据
    float *receiveData = new float[dataLength * nproc]; // 接收缓冲区

    // 调用自定义的 MPI_Allgather 函数
    mpi_allgatherSelfImpl(&sendData, dataLength, MPI_FLOAT, receiveData, dataLength, MPI_FLOAT, MPI_COMM_WORLD);

    // 打印所有接收到的数据以验证
    std::cout << "Process " << rank << " received data: ";
    for (int i = 0; i < dataLength * nproc; ++i) {
        std::cout << receiveData[i] << " ";
    }
    std::cout << std::endl;

    delete[] receiveData; // 释放动态分配的内存

    MPI_Finalize();
    return 0;
}
