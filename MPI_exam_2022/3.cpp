#include <iostream>
#include <mpi.h>

typedef struct {
    int m[3];
    float a[2];
    char c[5];
} mixtype;

void mpistruct(MPI_Datatype *newtype) {
    // 定义一个mixtype类型的变量，用于计算位移
    mixtype s;

    // 定义三个数组，分别存储数据类型、块长度和位移
    MPI_Datatype oldtype[3];
    int blocklen[3];
    MPI_Aint displaces[3];

    // 计算每个成员的地址并存储到displaces数组中
    MPI_Get_address(&s.m[0], &displaces[0]);
    MPI_Get_address(&s.a[0], &displaces[1]);
    MPI_Get_address(&s.c[0], &displaces[2]);

    // 计算相对位移
    displaces[1] -= displaces[0];
    displaces[2] -= displaces[0];
    displaces[0] = 0;

    // 设置每个块的长度
    blocklen[0] = 3;
    blocklen[1] = 2;
    blocklen[2] = 5;

    // 设置每个块的类型
    oldtype[0] = MPI_INT;
    oldtype[1] = MPI_FLOAT;
    oldtype[2] = MPI_CHAR;

    // 创建新的MPI数据类型
    MPI_Type_create_struct(3, blocklen, displaces, oldtype, newtype);
    MPI_Type_commit(newtype);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    // 定义并复制通信器
    MPI_Comm global_comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &global_comm);

    // 获取当前进程的秩和总进程数
    int myrank, numprocs;
    MPI_Comm_rank(global_comm, &myrank);
    MPI_Comm_size(global_comm, &numprocs);

    // 定义自定义数据类型
    MPI_Datatype newtp;
    mpistruct(&newtp);

    // 定义并初始化数据数组
    mixtype x[10];
    if (myrank == 0) {
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 3; ++j) x[i].m[j] = i + j;
            for (int j = 0; j < 2; ++j) x[i].a[j] = i * 0.1f + j;
            for (int j = 0; j < 5; ++j) x[i].c[j] = 'a' + j;
        }
    }

    // 进程0发送数据，进程1接收数据
    MPI_Status st;
    if (myrank == 0) {
        MPI_Send(x, 5, newtp, 1, 99, global_comm);
    } else if (myrank == 1) {
        MPI_Recv(x, 5, newtp, 0, 99, global_comm, &st);

        // 打印接收到的数据以验证
        for (int i = 0; i < 5; ++i) {
            std::cout << "Data element " << i << ":\n";
            std::cout << "m: ";
            for (int j = 0; j < 3; ++j) std::cout << x[i].m[j] << " ";
            std::cout << "\na: ";
            for (int j = 0; j < 2; ++j) std::cout << x[i].a[j] << " ";
            std::cout << "\nc: ";
            for (int j = 0; j < 5; ++j) std::cout << x[i].c[j];
            std::cout << "\n";
        }
    }

    // 释放自定义数据类型
    MPI_Type_free(&newtp);

    MPI_Finalize();
    return 0;
}
