#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

long long calculate(long long local_toss) {
  long long local_sum = 0;
  unsigned int seed = time(NULL);

  for (long long i = 0; i < local_toss; i++) {
    float x = (float)rand_r(&seed) / RAND_MAX;
    float y = (float)rand_r(&seed) / RAND_MAX;
    // in circle
    if (x * x + y * y <= 1)
      local_sum++;
  }
  return local_sum;
}

int main(int argc, char **argv) {
  // --- DON'T TOUCH ---
  MPI_Init(&argc, &argv);
  double start_time = MPI_Wtime();
  double pi_result;
  long long int tosses = atoi(argv[1]);
  int world_rank, world_size;
  // ---

  // TODO: init MPI
  // get size and rank
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  long long local_sum = 0, local_toss;
  if (world_rank > 0) {
    // TODO: handle workers
    local_toss = tosses / world_size;
    local_sum = calculate(local_toss);

    // send message
    MPI_Send(&local_sum, 1, MPI_LONG_LONG, 0, 0, MPI_COMM_WORLD);

  } else if (world_rank == 0) {
    // TODO: master
    local_toss = tosses - (world_size - 1) * (tosses / world_size);
    local_sum = calculate(local_toss);
  }

  if (world_rank == 0) {
    // TODO: process PI result
    long long tmp;
    MPI_Status status;
    for (int i = 1; i < world_size; i++) {
      MPI_Recv(&tmp, 1, MPI_LONG_LONG, i, 0, MPI_COMM_WORLD, &status);
      local_sum += tmp;
    }
    pi_result = 4 * local_sum / (double)tosses;

    // --- DON'T TOUCH ---
    double end_time = MPI_Wtime();
    printf("%lf\n", pi_result);
    printf("MPI running time: %lf Seconds\n", end_time - start_time);
    // ---
  }

  MPI_Finalize();
  return 0;
}
