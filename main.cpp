/* Copyright 2023 Stanford University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdio.h> // for printf
#include <stdlib.h> // for atoi
#include <cstring> // for memcpy
#include <algorithm> // for sort
#include <mpi.h>
#include <omp.h>

#define ROOT 3

// HEADERS
 #include <nccl.h>
// #include <rccl.h>
// #include <sycl.hpp>
// #include <ze_api.h>

// PORTS
 #define PORT_CUDA
// #define PORT_HIP
// #define PORT_SYCL

#include "comm.h"
namespace Comm = CommBench;

#include "coll.h"

// UTILITIES
#include "../CommBench/util.h"
void print_args();

// USER DEFINED TYPE
#define Type int
/*struct Type
{
  // int tag;
  int data[1];
  // complex<double> x, y, z;
};*/

int main(int argc, char *argv[])
{
  // INITIALIZE
  int myid;
  int numproc;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  int numthread;
  #pragma omp parallel
  if(omp_get_thread_num() == 0)
    numthread = omp_get_num_threads();
  // char machine_name[MPI_MAX_PROCESSOR_NAME];
  // int name_len = 0;
  // MPI_Get_processor_name(machine_name, &name_len);
  // printf("myid %d %s\n",myid, machine_name);

  // INPUT PARAMETERS
  if(argc != 6) {printf("arcgc: %d\n", argc); print_args(); MPI_Finalize(); return 0;}
  int library = atoi(argv[1]);
  int pattern = atoi(argv[2]);
  size_t count = atoi(argv[3]);
  int warmup = atoi(argv[4]);
  int numiter = atoi(argv[5]);

  // PRINT NUMBER OF PROCESSES AND THREADS
  if(myid == ROOT)
  {
    printf("\n");
    printf("Number of processes: %d\n", numproc);
    printf("Number of threads per proc: %d\n", numthread);
    printf("Number of warmup %d\n", warmup);
    printf("Number of iterations %d\n", numiter);

    printf("Pattern: %d\n", pattern);

    printf("Bytes per Type %lu\n", sizeof(Type));
    printf("Point-to-point (P2P) count %ld ( %ld Bytes)\n", count, count * sizeof(Type));
    printf("\n");
  }

  setup_gpu();

  // ALLOCATE
  Type *sendbuf = new Type[count * numproc];
  Type *recvbuf = new Type[count * numproc];
  Type *sendbuf_d;
  Type *recvbuf_d;
#ifdef PORT_CUDA
  cudaMalloc(&sendbuf_d, count * sizeof(Type) * numproc);
  cudaMalloc(&recvbuf_d, count * sizeof(Type) * numproc);
#elif defined PORT_HIP
  hipMalloc(&sendbuf_d, count * sizeof(Type) * numproc);
  hipMalloc(&recvbuf_d, count * sizeof(Type) * numproc);
#elif defined PORT_SYCL
  sycl::queue q(sycl::gpu_selector_v);
  sendbuf_d = sycl::malloc_device<Type>(count * numproc, q);
  recvbuf_d = sycl::malloc_device<Type>(count * numproc, q);
#endif

  // Coll::Comm_direct<Type> coll(MPI_COMM_WORLD, Comm::NCCL);
  // Coll::Comm_direct<Type> coll(MPI_COMM_WORLD, 4, Comm::NCCL, Comm::IPC);
  // Coll::Comm_direct<Type> coll(MPI_COMM_WORLD, 4, 1, Comm::NCCL, Comm::IPC, Comm::IPC);
  Comm::Comm<Type> coll(MPI_COMM_WORLD, (Comm::library) library);

  switch(pattern) {
    case 0:
      if(myid == ROOT)
        printf("TEST P2P\n");
      coll.add(sendbuf_d, 0, recvbuf_d, 0, count, 0, 1);
      break;
    case 1:
      if(myid == ROOT)
        printf("TEST GATHER\n");
      for(int p = 0; p < numproc; p++)
        coll.add(sendbuf_d, 0, recvbuf_d, p * count, count, p, ROOT);
      break;
    case 2:
      if(myid == ROOT)
        printf("TEST SCATTER\n");
      for(int p = 0; p < numproc; p++)
        coll.add(sendbuf_d, p * count, recvbuf_d, 0, count, ROOT, p);
      break;
    case 3:
      if(myid == ROOT)
        printf("TEST REDUCE\n");
      break;
    case 4:
      if(myid == ROOT)
        printf("TEST BROADCAST\n");
      for(int p = 0; p < numproc; p++)
        coll.add(sendbuf_d, 0, recvbuf_d, 0, count, ROOT, p);
      break;
    case 5:
      if(myid == ROOT)
        printf("TEST ALL-TO-ALL\n");
      for(int sender = 0; sender < numproc; sender++)
        for(int recver = 0; recver < numproc; recver++)
          coll.add(sendbuf_d, recver * count, recvbuf_d, sender * count, count, sender, recver);
      break;
    case 6:
      if(myid == ROOT)
        printf("TEST ALL-REDUCE\n");
      break;
    case 7:
      if(myid == ROOT)
        printf("TEST ALL-GATHER\n");
      for(int sender = 0; sender < numproc; sender++)
        for(int recver = 0; recver < numproc; recver++)
          coll.add(sendbuf_d, 0, recvbuf_d, sender * count, count, sender, recver);
      break;
    case 8:
      if(myid == ROOT)
        printf("TEST REDUCE-SCATTER\n");
      break;
  }
  // coll.init();


  for(int iter = 0; iter < numiter; iter++)
    validate(sendbuf_d, recvbuf_d, count, pattern, coll);

  coll.measure(warmup, numiter);

  return 0;


// DEALLOCATE
#ifdef PORT_CUDA
  cudaFree(sendbuf_d);
  cudaFree(recvbuf_d);
#elif defined PORT_HIP
  hipFree(sendbuf_d);
  hipFree(recvbuf_d);
#elif defined PORT_SYCL
  sycl::free(sendbuf_d, q);
  sycl::free(recvbuf_d, q);
#else
  delete[] sendbuf_d;
  delete[] recvbuf_d;
#endif

  // FINALIZE
  MPI_Finalize();

  return 0;
} // main()

void print_args() {

  int myid;
  int numproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);

  if(myid == ROOT) {
    printf("\n");
    printf("CollBench requires four arguments:\n");
    printf("1. library:\n");
    printf("      0 for IPC\n");
    printf("      1 for MPI\n");
    printf("      2 for NCCL\n");
    printf("2. pattern:\n");
    printf("      1 for Gather\n");
    printf("      2 for Scatter\n");
    printf("      3 for Reduce\n");
    printf("      4 for Broadcast\n");
    printf("      5 for Alltoall\n");
    printf("      6 for Allreduce\n");
    printf("      7 for Allgather\n");
    printf("      8 for ReduceScatter\n");
    printf("3. count: number of 4-byte elements\n");
    printf("4. warmup: number of warmup rounds\n");
    printf("5. numiter: number of measurement rounds\n");
    printf("where on can run CollBench as\n");
    printf("mpirun ./CollBench count warmup numiter\n");
    printf("\n");
  }
}

