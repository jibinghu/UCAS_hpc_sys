#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    MPI_Datatype indexed_type;
    int block_lengths[2] = {2, 3};
    int block_displacements[2] = {1, 5};  // 基于 MPI_INT 假定是 sizeof(int)
    MPI_Type_indexed(2, block_lengths, block_displacements, MPI_INT, &indexed_type);
    MPI_Type_commit(&indexed_type);

    MPI_Aint extent;
    MPI_Type_extent(indexed_type, &extent);
    // int MPI_Type_get_extent(MPI_Datatype, MPI_Aint*, MPI_Aint*)
    printf("Extent of indexed_type: %ld bytes.\n", extent);

    MPI_Type_free(&indexed_type);
    MPI_Finalize();
    return 0;
}
