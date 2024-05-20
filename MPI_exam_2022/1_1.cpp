#include <mpi.h>
#include <iostream>


void create_2d_communicators(int p, int q, MPI_Comm &rowComm, MPI_Comm &colComm){

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    
    if (p * q != size) {
        if (rank == 0){
        std::cerr << "Correct your program" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    int rowColor = rank / q;
    // int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
    MPI_Comm_split(MPI_COMM_WORLD, rowColor, rank, &rowComm);
    int colColor = rank % q;
    MPI_Comm_split(MPI_COMM_WORLD, colColor, rank, &colComm);
}

int main(int argc, char** argv) {
    // 初始化 MPI 环境
    MPI_Init(NULL, NULL);

    MPI_Comm rowComm,colComm;
    create_2d_communicators(3,3,rowComm,colComm);

    int rowRank, colRank, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_rank(rowComm, &rowRank);
    MPI_Comm_rank(colComm, &colRank);

    std::cout << "rank of world: " << rank << std::endl; 
    std::cout << "rank of row communicator: " << rowRank << std::endl;
    std::cout << "ranl of col communicator: " << colRank << std::endl;

    MPI_Comm_free(&rowComm);
    MPI_Comm_free(&colComm);
    // 释放 MPI 的一些资源
    MPI_Finalize();
}