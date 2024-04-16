#include <string>
#include "mpi.h"
#include <chrono>
#include <cmath>

using namespace std;

void printArray(int* arr, int size) {
    for (int i = 0; i < size; ++i) {
        fprintf(stdout, "%d ", arr[i]);
    }
    fprintf(stdout, "\n");
    fflush(stdout);
}

int* oddEvenSort(int* arr, int N, int rank, int size){
//    fprintf(stdout, "Start\n");
//    fflush(stdout);
    double startTime;
    if (rank == 0)
    {
        startTime = MPI_Wtime();
    }
    for (int i = 0; i < N+1; i++) {
        if (i % 2 == 0) {
            if (rank != 0) {
                for (int j = (rank - 1) * 2; j < N-1; j += (size - 1) * 2) {
                    if (arr[j] > arr[j + 1]) {
                        int temp = arr[j + 1];
                        arr[j + 1] = arr[j];
                        arr[j] = temp;
                    }
                    MPI_Send(&arr[j], 1, MPI_INT, 0,j+1, MPI_COMM_WORLD);
                    MPI_Send(&arr[j + 1], 1, MPI_INT, 0,j+2, MPI_COMM_WORLD);
                }
            } else {
                int len = N;
                if (N%2 == 1) len--;
                for (int i = 0; i < len; i++){
                    MPI_Status status;
                    MPI_Recv(&arr[i], 1, MPI_INT, MPI_ANY_SOURCE, i+1, MPI_COMM_WORLD, &status);
                }
            }
        } else {
            if (rank != 0) {
                for (int j = (rank - 1) * 2 + 1; j < N-1; j += (size - 1) * 2) {
                    if (arr[j] > arr[j + 1]) {
                        int temp = arr[j + 1];
                        arr[j + 1] = arr[j];
                        arr[j] = temp;
                    }
                    MPI_Send(&arr[j], 1, MPI_INT, 0,j+1, MPI_COMM_WORLD);
                    MPI_Send(&arr[j + 1], 1, MPI_INT, 0,j+2, MPI_COMM_WORLD);
                }
            } else {
                int len = N;
                if (N%2 == 0) len--;
                for (int i = 1; i < len; i++){
                    MPI_Status status;
                    MPI_Recv(&arr[i], 1, MPI_INT, MPI_ANY_SOURCE, i+1, MPI_COMM_WORLD, &status);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
//            printArray(arr,N);
//            fprintf(stdout, "%d\n", i);
//            fflush(stdout);
            for (int i = 1; i < size; i++){
                MPI_Send(arr, N, MPI_INT, i,0, MPI_COMM_WORLD);
            }
        }
//        MPI_Barrier(MPI_COMM_WORLD);
        if (rank != 0) {
//            fprintf(stdout, "Process %d before\n", rank);
//            printArray(arr,N);
//            fflush(stdout);
            MPI_Status status;
            MPI_Recv(arr, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
//            fprintf(stdout, "Process %d after\n", rank);
//            printArray(arr,N);
//            fflush(stdout);
        }
    }
    if (rank == 0)
    {
        double endTime = MPI_Wtime();
        double elapsedTime = endTime - startTime;
        fprintf(stdout, "Elapsed parallel time:  %f\n", elapsedTime);
        fflush(stdout);
    }
    return arr;
}

int* oddEvenSortNonParallel(int* arr, int N){
    double startTime;
    startTime = MPI_Wtime();
    for (int i = 0; i < N; i++) {
        for (int pos = 0; pos < N; pos++) {
            if (i % 2 == 0) {
                if (pos % 2 == 0) {
                    if (arr[pos] > arr[pos + 1]) {
                        int temp = arr[pos + 1];
                        arr[pos + 1] = arr[pos];
                        arr[pos] = temp;
                    }
                }
            } else {
                if (pos % 2 == 1) {
                    if (arr[pos] > arr[pos + 1]) {
                        int temp = arr[pos + 1];
                        arr[pos + 1] = arr[pos];
                        arr[pos] = temp;
                    }
                }
            }
        }
    }
    double endTime = MPI_Wtime();
    double elapsedTime = endTime - startTime;
    fprintf(stdout, "Elapsed nonparallel time:  %f\n", elapsedTime);
    fflush(stdout);
    return arr;
}


double f(double a)
{
    return (4.0 / (1.0+ a*a));
}

int main(int argc, char **argv) {
    int N;
    N = 2000;
    srand(time(nullptr));
    int myArray[N];

    for (int i = 0; i < N; ++i) {
        myArray[i] = 0 + rand() % (10000 - 0 + 1);
    };
    int rank, size;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(processor_name,&namelen);

    int* sortedArray = oddEvenSort(myArray,N, rank, size);
    if (rank == 0) {
//        printArray(sortedArray, N);
    }

    if (rank == 0) {
        sortedArray = oddEvenSortNonParallel(myArray, N);
//        printArray(sortedArray, N);
    }
    MPI_Finalize();

    return 0;

}
