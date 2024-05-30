#include <stdio.h>
#include "mpi.h"

#define row_P 4 // 行处理器数
#define col_P 3 // 列处理器数

#define row_A 9 // 矩阵A的行数
#define col_A 7 // 矩阵A的列数

// 用卷帘方式实现8*6的矩阵放至4*3的处理器上
// 注意这里的处理器分配矩阵不是传统的最后一个处理器多余处理，而是前i个处理器每个处理器负责1行/列

void rowcolcomm(int myid, MPI_Comm comm)
{
    int rowid, colid;      
    int ma, ka;            
    int rowcolor, colcolor; 
    int i, j;
    int A[row_A][col_A];   

    MPI_Comm rowcomm, colcomm;

    // 计算行颜色和列颜色，用于分割通信器
    rowcolor = myid / col_P;
    MPI_Comm_split(comm, rowcolor, myid, &rowcomm); 
    MPI_Comm_rank(rowcomm, &colid);

    colcolor = myid % col_P;
    MPI_Comm_split(comm, colcolor, myid, &colcomm); 
    MPI_Comm_rank(colcomm, &rowid);

    // 计算每个进程负责的矩阵块的行数
    if (rowid < row_A % row_P) {
        ma = row_A / row_P + 1;
    } else {
        ma = row_A / row_P;
    }

    // 计算每个进程负责的矩阵块的列数
    if (colid < col_A % col_P) {
        ka = col_A / col_P + 1;
    } else {
        ka = col_A / col_P;
    }

    // 打印进程ID和负责的矩阵块大小
    printf("Process %3d ma=%d ka=%d Aij is ", myid, ma, ka);
    for (i = 0; i < ma; i++)
    {
        for (j = 0; j < ka; j++)
        {
            // 计算并存储矩阵块的元素
            A[i][j] = i * row_P + j * col_P + rowid + colid;
            printf("%3d ", A[i][j]);
        }
    }
    MPI_Barrier(rowcomm);
    MPI_Barrier(colcomm);
    printf("\n");
}

int main(int argc, char *argv[])
{
    int rank; // 当前进程的rank
    int size; // 总进程数
    MPI_Comm mycomm;

    // 初始化MPI环境
    MPI_Init(&argc, &argv);
    MPI_Comm_dup(MPI_COMM_WORLD, &mycomm);
    MPI_Comm_rank(mycomm, &rank);
    MPI_Comm_size(mycomm, &size);

    // 检查总进程数是否匹配
    if (size != row_P * col_P) {
        if (rank == 0) {
            printf("Error: Number of processes must be %d\n", row_P * col_P);
        }
        MPI_Finalize();
        return -1;
    }

    // 调用分割行和列通信器的函数
    rowcolcomm(rank, mycomm);

    // 结束MPI环境
    MPI_Finalize();
    return 0;
}
