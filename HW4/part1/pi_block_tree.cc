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

  // TODO: MPI init
  // get size and rank
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // TODO: binary tree redunction
  long long local_sum = 0, local_toss;
  if (world_rank == 0)
    local_toss = tosses - (world_size - 1) * (tosses / world_size);
  else
    local_toss = tosses / world_size;

  local_sum = calculate(local_toss);

  for (int k = 2; k <= world_size; k *= 2) {
    if (world_rank % k != 0) { // sender
      MPI_Send(&local_sum, 1, MPI_LONG_LONG, world_rank - k / 2, 0,
               MPI_COMM_WORLD);
      break;
    } else { // receicer
      long long rec_count;
      MPI_Status status;
      MPI_Recv(&rec_count, 1, MPI_LONG_LONG, world_rank + k / 2, 0,
               MPI_COMM_WORLD, &status);
      local_sum += rec_count;
    }
  }

  if (world_rank == 0) {
    // TODO: PI result
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
