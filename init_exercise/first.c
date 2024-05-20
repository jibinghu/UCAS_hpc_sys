#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // 动态分配数组A的内存
    int* A = (int*)malloc(4 * sizeof(int)); // 假设数组长度为4

    // 初始化数组
    if(world_rank == 0) {
        for(int i = 0; i < 4; i++) A[i] = 1; // 为进程0的数组元素赋值为1
    } else if (world_rank == 1) {
        for(int i = 0; i < 4; i++) A[i] = 2; // 为进程1的数组元素赋值为2
    }

    // 计算每个进程的局部总和
    int local_sum = 0;
    for(int i = 0; i < 4; i++) {
        local_sum += A[i];
    }

    if(world_rank == 1) {
        // 进程1向进程0发送局部总和
        MPI_Send(&local_sum, 1, MPI_INT, 0, 99, MPI_COMM_WORLD);
    }

    if(world_rank == 0) {
        int received_sum;
        MPI_Status status;
        // 接收进程1的局部总和
        MPI_Recv(&received_sum, 1, MPI_INT, 1, 99, MPI_COMM_WORLD, &status);
        local_sum += received_sum; // 汇总总和
        printf("Total sum: %d\n", local_sum); // 输出总和
    }

    // 释放动态分配的内存
    free(A);

    MPI_Finalize();
    return 0;
}
