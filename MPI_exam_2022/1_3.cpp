#include <mpi.h>
#include <iostream>


// 在1_1的基础上进行1_3的code
void create_2d_communicators(int p, int q, MPI_Comm &rowComm, MPI_Comm &colComm){

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    if (p * q != size) {
        if (rank == 0){
        std::cerr << "Correct your program" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    int rowColor = rank / q;
    // int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
    MPI_Comm_split(MPI_COMM_WORLD, rowColor, rank, &rowComm);
    int colColor = rank % q;
    MPI_Comm_split(MPI_COMM_WORLD, colColor, rank, &colComm);
}

void bcast_from_p00(int &p_00, MPI_Comm rowComm, MPI_Comm colComm){

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int rowRank;
    MPI_Comm_rank(rowComm, &rowRank);
    int colRank;
    MPI_Comm_rank(colComm, &colRank);

    // 先将P_00广播至每行的第一列
    // if (rowRank == 0){
       if (world_rank == 0){
        // int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
        MPI_Bcast(&p_00, 1, MPI_INT, 0, colComm);
    }

    // 再将每行第一列中的p_00广播至所有列
    MPI_Bcast(&p_00, 1, MPI_INT, 0, rowComm);

}

int main(int argc, char** argv) {
    // 初始化 MPI 环境
    MPI_Init(NULL, NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm rowComm,colComm;
    create_2d_communicators(3,3,rowComm,colComm);

    int p_00;
    if (rank == 0)
        p_00 = 520; // 今天是2024年5月20日，彩蛋：祝天下有情人终成眷属！

    bcast_from_p00(p_00, rowComm, colComm);

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "My process ID: " << rank << ", and my p_00 value is: " << p_00 << std::endl;

    // int rowRank, colRank, rank;
    // MPI_Comm_rank(rowComm, &rowRank);
    // MPI_Comm_rank(colComm, &colRank);

    // std::cout << "rank of world: " << rank << std::endl; 
    // std::cout << "rank of row communicator: " << rowRank << std::endl;
    // std::cout << "ranl of col communicator: " << colRank << std::endl;

    MPI_Comm_free(&rowComm);
    MPI_Comm_free(&colComm);
    // 释放 MPI 的一些资源
    MPI_Finalize();
}