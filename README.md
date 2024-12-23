# Parallel-Programming-Final

> *Optimizing the Traveling Salesman Problem with Parallel Algorithms*\
> Group Member: 蔡承翰、蔡師睿、林書楷

## Run GA or ACO with Non-parallel, OpenMP, Thread, or MPI Implementations
```bash
$ make run testcase=<path_to_testcase_file> algo=<algorithm_type> parallel=<parallel_type>

# Replace <path_to_testcase_file> with the path to the testcase file, ex: tsp_graph/a280.tsp
# Replace <algorithm_type> with the algorithm type (aco or ga)
# Replace <parallel_type> with the parallel type (serial, omp, thread, or mpi)
#
# Example:
# make run testcase=./tsp_graph/a280.tsp algo=ga parallel=thread
```

## Run ACO with CUDA
```bash
$ make cuda testcase=<path_to_testcase_file>

# Replace <path_to_testcase_file> with the path to the testcase file, ex: tsp_graph/a280.tsp
#
# Example:
# make cuda testcase=./tsp_graph/a280.tsp
```

> [!NOTE]
> Due to some unknown issue, the cuda code cannot run on the hpc1 but can run on my server.
