#include <mpi.h>
#include <iostream>
using namespace std;

void make_newtype(int m, int n, int lda, MPI_Datatype *newtype) {
    MPI_Datatype vec_type;
    const int count = m;
    const int length = n;
    const int stride = lda;

    MPI_Type_vector(count, length, stride - m, MPI_INT, &vec_type);
    MPI_Type_create_resized(vec_type, 0, sizeof(int) * lda * m * 2, newtype);
    MPI_Type_commit(newtype);
    MPI_Type_free(&vec_type);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int m = 2; // 块矩阵的行数
    int n = 2; // 块矩阵的列数
    int lda = 6; // 原始矩阵的列数
    MPI_Datatype newtype;
    make_newtype(m, n, lda, &newtype);

    if (rank == 0) {
        int A[6][6] = {
            {1, 2, 3, 4, 0, 0},
            {5, 6, 7, 8, 0, 0},
            {9, 10, 11, 12, 0 ,0},
            {13, 14, 15, 16, 0, 0},
            {17, 18, 19, 20, 0, 0},
            {21, 22, 23, 24, 0, 0}
        };

        // 发送 A00 和 A20 块矩阵
        MPI_Send(&A[0][0], 2, newtype, 1, 0, MPI_COMM_WORLD);
        // MPI_Send(&A[4][0], 1, newtype, 1, 1, MPI_COMM_WORLD);
    } else if (rank == 1) {
        int A_recv[2][2];
        int B_recv[2][2];

        // 接收 A00 和 A20 块矩阵
        MPI_Recv(&A_recv[0][0], 2, newtype, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // MPI_Recv(&B_recv[0][0], 1, newtype, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        cout << "Received A_00:" << endl;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                cout << A_recv[i][j] << " ";
            }
            cout << endl;
        }

        cout << "Received A_20:" << endl;
        for (int i = 2; i < 2+m; i++) {
            for (int j = 2; j < 2+n; j++) {
                cout << B_recv[i][j] << " ";
            }
            cout << endl;
        }
    }

    MPI_Type_free(&newtype);
    MPI_Finalize();

    return 0;
}
