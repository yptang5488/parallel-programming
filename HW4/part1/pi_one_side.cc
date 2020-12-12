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

  MPI_Win win;

  // TODO: MPI init
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  long long *total_sum;
  long long local_sum;

  if (world_rank == 0) {
    // Master
    // initialize
    MPI_Alloc_mem(sizeof(long long), MPI_INFO_NULL, &total_sum);
    *total_sum = 0;
    MPI_Win_create(total_sum, sizeof(long long), sizeof(long long),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    local_sum = calculate(tosses / world_size +
                          (tosses % world_size <= world_rank ? 0 : 1));
  } else {
    // Workers
    MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    local_sum = calculate(tosses / world_size +
                          (tosses % world_size <= world_rank ? 0 : 1));

    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
    MPI_Accumulate(&local_sum, 1, MPI_LONG_LONG, 0, 0, 1, MPI_LONG_LONG,
                   MPI_SUM, win);
    MPI_Win_unlock(0, win);
  }

  MPI_Win_fence(0, win);
  MPI_Win_free(&win);

  if (world_rank == 0) {
    // TODO: handle PI result
    local_sum += *total_sum;
    pi_result = 4 * local_sum / (double)tosses;

    MPI_Free_mem(total_sum);

    // --- DON'T TOUCH ---
    double end_time = MPI_Wtime();
    printf("%lf\n", pi_result);
    printf("MPI running time: %lf Seconds\n", end_time - start_time);
    // ---
  }

  MPI_Finalize();
  return 0;
}