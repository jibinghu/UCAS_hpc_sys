#include <mpi.h>
#include <stdio.h>

typedef struct {
    double weight;
    int count;
    char label;
} Item;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    MPI_Datatype item_type;
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_INT, MPI_CHAR};
    int block_lengths[3] = {1, 1, 1};
    MPI_Aint offsets[3];

    MPI_Get_address(&((Item*)0)->weight, &offsets[0]);
    MPI_Get_address(&((Item*)0)->count, &offsets[1]);
    MPI_Get_address(&((Item*)0)->label, &offsets[2]);

    MPI_Type_create_struct(3, block_lengths, offsets, types, &item_type);
    MPI_Type_commit(&item_type);

    MPI_Aint extent;
    MPI_Type_extent(item_type, &extent);
    // int MPI_Type_get_extent(MPI_Datatype, MPI_Aint*, MPI_Aint*)
    printf("Extent of item_type: %ld bytes.\n", extent);

    MPI_Type_free(&item_type);
    MPI_Finalize();
    return 0;
}
