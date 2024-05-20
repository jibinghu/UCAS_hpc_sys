#include <mpi.h>
#include <stdio.h>
#include <cstddef>

typedef struct {
    double weight;
    int count;
} Item;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    // 创建 MPI 类型
    MPI_Datatype item_type;
    MPI_Datatype type[2] = {MPI_DOUBLE, MPI_INT};
    int blocklen[2] = {1, 1};
    MPI_Aint disp[2];

    disp[0] = offsetof(Item, weight);
    disp[1] = offsetof(Item, count);

    MPI_Type_create_struct(2, blocklen, disp, type, &item_type);
    MPI_Type_commit(&item_type);

    // 获取类型的扩展
    MPI_Aint extent;
    MPI_Type_extent(item_type, &extent);
    printf("The extent of the item_type is %ld bytes.\n", extent);

    MPI_Type_free(&item_type);
    MPI_Finalize();
    return 0;
}
