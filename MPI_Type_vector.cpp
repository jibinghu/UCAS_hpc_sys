#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    MPI_Datatype vector_type;
    int block_length = 2;
    int stride = 3;
    int count = 4;
    MPI_Type_vector(count, block_length, stride, MPI_DOUBLE, &vector_type);
    MPI_Type_commit(&vector_type);

    MPI_Aint extent;
    MPI_Type_extent(vector_type, &extent);
    // int MPI_Type_get_extent(MPI_Datatype, MPI_Aint*, MPI_Aint*)
    printf("Extent of vector_type: %ld bytes.\n", extent);

    MPI_Type_free(&vector_type);
    MPI_Finalize();
    return 0;
}
