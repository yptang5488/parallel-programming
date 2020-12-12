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
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  long long total_sum = 0, local;
  local = tosses / world_size + (tosses % world_size <= world_rank ? 0 : 1);
  total_sum = calculate(local);

  // TODO: use MPI_Reduce
  MPI_Reduce(&local, &total_sum, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if (world_rank == 0) {
    // TODO: PI result
    pi_result = 4 * total_sum / (double)tosses;

    // --- DON'T TOUCH ---
    double end_time = MPI_Wtime();
    printf("%lf\n", pi_result);
    printf("MPI running time: %lf Seconds\n", end_time - start_time);
    // ---
  }

  MPI_Finalize();
  return 0;
}
