#include <vector>

namespace Coll {

  template <typename T>
  struct P2P_h {
    public:
    T *sendbuf;
    const size_t sendoffset;
    T *recvbuf;
    const size_t recvoffset;
    const size_t count;
    const int sendid;
    const int recvid;

    P2P_h(T *sendbuf, size_t sendoffset, T *recvbuf, size_t recvoffset, size_t count, int sendid, int recvid)
    : sendbuf(sendbuf), sendoffset(sendoffset), recvbuf(recvbuf), recvoffset(recvoffset), count(count), sendid(sendid), recvid(recvid) {};

  };

  template <typename T>
  class Comm_h {

    protected:

    const MPI_Comm &comm_mpi;
    int groupsize;
    int subgroupsize;
    Comm::library selflib;
    Comm::library intralib;
    Comm::library interlib;

    std::vector<P2P_h<T>*> level_2;
    std::vector<P2P_h<T>*> level_1;
    std::vector<P2P_h<T>*> level_0;

    public:

    Comm_h(const MPI_Comm &comm_mpi, Comm::library lib)
    : Comm_h(comm_mpi, 0, 0, lib, lib, lib) {};
    Comm_h(const MPI_Comm &comm_mpi, int groupsize, Comm::library inter, Comm::library intra)
    : Comm_h(comm_mpi, groupsize, 0, inter, intra, intra) {};
    Comm_h(const MPI_Comm &comm_mpi, int groupsize, int subgroupsize, Comm::library interlib, Comm::library intralib, Comm::library selflib)
    : comm_mpi(comm_mpi), groupsize(groupsize), subgroupsize(subgroupsize), interlib(interlib), intralib(intralib), selflib(selflib) {}

    void add(T *sendbuf, size_t sendoffset, T *recvbuf, size_t recvoffset, size_t count, int sendid, int recvid) {

      int myid;
      int numproc;
      MPI_Comm_rank(comm_mpi, &myid);
      MPI_Comm_size(comm_mpi, &numproc);

      P2P_h<T> *p2p_temp = new P2P_h<T>(sendbuf, sendoffset, recvbuf, recvoffset, count, sendid, recvid);
      if((subgroupsize > 0) && (sendid / subgroupsize == recvid / subgroupsize))
        level_2.push_back(p2p_temp);
      else if((groupsize > 0) && (sendid / groupsize == recvid / groupsize))
        level_1.push_back(p2p_temp);
      else
        level_0.push_back(p2p_temp);
    };

    ~Comm_h() {
      for(P2P_h<T> *p2p : this->level_2)
        delete p2p;
      for(P2P_h<T> *p2p : this->level_1)
        delete p2p;
      for(P2P_h<T> *p2p : this->level_0)
        delete p2p;
    }
  };

  template<typename T>
  class Comm_direct : public Comm_h<T> {

    Comm::Comm<T> *self = new Comm::Comm<T>(this->comm_mpi, this->selflib);
    Comm::Comm<T> *intra = new Comm::Comm<T>(this->comm_mpi, this->intralib);
    Comm::Comm<T> *inter = new Comm::Comm<T>(this->comm_mpi, this->interlib);

    public:

    using Comm_h<T>::Comm_h;

    void init() {

      int myid;
      int numproc;
      MPI_Comm_rank(this->comm_mpi, &myid);
      MPI_Comm_size(this->comm_mpi, &numproc);

      if(myid == ROOT)
        printf("self\n");
      for(P2P_h<T> *p2p : this->level_2)
        self->add(p2p->sendbuf, p2p->sendoffset, p2p->recvbuf, p2p->recvoffset, p2p->count, p2p->sendid, p2p->recvid);

      if(myid == ROOT)
        printf("intra\n");
      for(P2P_h<T> *p2p : this->level_1)
        intra->add(p2p->sendbuf, p2p->sendoffset, p2p->recvbuf, p2p->recvoffset, p2p->count, p2p->sendid, p2p->recvid);

      if(myid == ROOT)
        printf("inter\n");
      for(P2P_h<T> *p2p : this->level_0)
        inter->add(p2p->sendbuf, p2p->sendoffset, p2p->recvbuf, p2p->recvoffset, p2p->count, p2p->sendid, p2p->recvid);
    };
    ~Comm_direct() {
      delete self;
      delete intra;
      delete inter;
    }
    void measure(int warmup, int numiter) {
      self->measure(warmup, numiter);
      intra->measure(warmup, numiter);
      inter->measure(warmup, numiter);
    }
    void run() {
      inter->start();
      intra->start();
      self->start();
      self->wait();
      intra->wait();
      inter->wait();
    }
  };

  template<typename T>
  class Comm_striped : public Comm_h<T> {

    Comm::Comm<T> *self = new Comm::Comm<T>(this->comm_mpi, this->selflib);
    Comm::Comm<T> *intra = new Comm::Comm<T>(this->comm_mpi, this->intralib);
    Comm::Comm<T> *partition = new Comm::Comm<T>(this->comm_mpi, this->intralib);
    Comm::Comm<T> *translation = new Comm::Comm<T>(this->comm_mpi, this->interlib);
    Comm::Comm<T> *assembly = new Comm::Comm<T>(this->comm_mpi, this->intralib);

    std::vector<size_t> sendcount;
    std::vector<size_t> recvcount;
    std::vector<T*> sendbuf;
    std::vector<T*> recvbuf;

    public:

    Comm_striped(const MPI_Comm &comm_mpi, int groupsize, int subgroupsize, Comm::library interlib, Comm::library intralib, Comm::library selflib) 
    : Comm_h<T>::Comm_h(comm_mpi, groupsize, subgroupsize, interlib, intralib, selflib) {};

    void init() {

      int myid;
      int numproc;
      MPI_Comm_rank(this->comm_mpi, &myid);
      MPI_Comm_size(this->comm_mpi, &numproc);

      if(myid == ROOT)
        printf("self\n");
      for(P2P_h<T> *p2p : this->level_2)
        self->add(p2p->sendbuf, p2p->sendoffset, p2p->recvbuf, p2p->recvoffset, p2p->count, p2p->sendid, p2p->recvid);

      if(myid == ROOT)
        printf("intra\n");
      for(P2P_h<T> *p2p : this->level_1)
        intra->add(p2p->sendbuf, p2p->sendoffset, p2p->recvbuf, p2p->recvoffset, p2p->count, p2p->sendid, p2p->recvid);

      if(myid == ROOT)
        printf("partition\n");
      {
        for(P2P_h<T> *p2p : this->level_0) {
          int sendgroup = p2p->sendid / this->groupsize;
          int recvgroup = p2p->recvid / this->groupsize;
          int partial[this->groupsize];
          for(int p = 0; p < this->groupsize; p++) {
            partial[p] = p2p->count / this->groupsize;
          }
	}
            
       /*     for(int recver = (sendid / this->groupsize) * groupsize; p < (sendid / this->group 
            sendbuf.push_back(new T*[count / this->groupsize]);
	    partition->add(p2p->sendbuf, myid % this->groupsize * count / this->groupsize, sendbuf.back, 0, count / this->groupsize, sendid, a
	  }

        for(P2P_h<T> *p2p : this->level_0)
          if(myid / this->groupsize == p2p->recvid / this->groupsize)
            recvcount.push_back(p2p->count / this->groupsize);

	for(int count : sendcount)

	for(int count : recvcount)
          recvbuf.push_back(new T*[count]);


	printf("myid %d numsend_temp %d numrecv_temp %d\n", myid, numsend, numrecv); */
      }

      if(myid == ROOT)
        printf("translation\n");
      for(P2P_h<T> *p2p : this->level_0)
        translation->add(p2p->sendbuf, p2p->sendoffset, p2p->recvbuf, p2p->recvoffset, p2p->count, p2p->sendid, p2p->recvid);

      if(myid == ROOT)
        printf("assembly\n");
    };
    void measure(int warmup, int numiter) {
      self->measure(warmup, numiter);
      intra->measure(warmup, numiter);
      partition->measure(warmup, numiter);
      translation->measure(warmup, numiter);
      assembly->measure(warmup, numiter);
    }
    void run() {
      partition->start();
      partition->wait();
      translation->start();
      intra->start();
      self->start();
      self->wait();
      intra->wait();
      translation->wait();
      assembly->start();
      assembly->wait();
    }
  };

}


template <typename T>
class Comm_striped {

  Comm::Comm<T> *self;
  Comm::Comm<T> *intra;
  Comm::Comm<T> *partition;
  Comm::Comm<T> *translation;
  Comm::Comm<T> *assembly;

  int numsend = 0;
  int numrecv = 0;
  T **sendbuf;
  T **recvbuf;
  MPI_Comm comm_mpi;
  int groupsize;
  int subgroupsize;

  public:

  Comm_striped(const MPI_Comm &comm_mpi, int groupsize, int subgroupsize, Comm::library interlib, Comm::library intralib, Comm::library selflib)
    : comm_mpi(comm_mpi), groupsize(groupsize), subgroupsize(subgroupsize) {
    int myid;
    int numproc;
    MPI_Comm_rank(comm_mpi, &myid);
    MPI_Comm_size(comm_mpi, &numproc);

    if(myid == ROOT)
      printf("Create Comm_striped with p = %d, g = %d, k = %d\n", numproc, groupsize, subgroupsize);

    if(myid == ROOT) printf("level 2 ");
    self = new Comm::Comm<T>(comm_mpi, selflib);
    if(myid == ROOT) printf("level 1 ");
    intra = new Comm::Comm<T>(comm_mpi, intralib);
    if(myid == ROOT) printf("partition ");
    partition = new Comm::Comm<T>(comm_mpi, intralib);
    if(myid == ROOT) printf("translation ");
    translation = new Comm::Comm<T>(comm_mpi, interlib);
    if(myid == ROOT) printf("assembly ");
    assembly = new Comm::Comm<T>(comm_mpi, intralib);
  }

  void add(T *sendbuf, size_t sendoffset, T *recvbuf, size_t recvoffset, size_t count, int sendid, int recvid) {
    int myid;
    int numproc;
    MPI_Comm_rank(comm_mpi, &myid);
    MPI_Comm_size(comm_mpi, &numproc);

    if(sendid / subgroupsize == recvid / subgroupsize) {
      if(myid == ROOT) printf("level 2 ");
      self->add(sendbuf, sendoffset, recvbuf, recvoffset, count, sendid, recvid);
    }
    else {
      if(sendid / groupsize == recvid / groupsize) {
        if(myid == ROOT) printf("level 1 ");
        intra->add(sendbuf, sendoffset, recvbuf, recvoffset, count, sendid, recvid);
      }
      else {
        if(myid == ROOT) printf("level 0 ");
        if(myid / groupsize == sendid / groupsize) {
          T **sendbuf_temp = new T*[numsend + 1];
          memcpy(sendbuf_temp, sendbuf, numsend * sizeof(T*));
          delete[] sendbuf_temp;
          sendbuf = sendbuf_temp;
          if(myid == sendid)
            sendbuf[numsend + 1] = new T[count];
          numsend++;
        }
        if(myid == ROOT) printf("partition ");
        int sendgroup = sendid / groupsize;
        for(int p = 0; p < groupsize; p++)
          for(int recver = sendgroup * groupsize; recver < (sendgroup + 1) * groupsize; recver++)
            if(recver != sendid)
              partition.add(sendbuf, sendoffset,        count, sendid, recver);
        translation->add(sendbuf, sendoffset, recvbuf, recvoffset, count, sendid, recvid);
      }
    }
  }
  void measure(int warmup, int numiter) {
    self->measure(warmup, numiter);
    intra->measure(warmup, numiter);
    translation->measure(warmup, numiter);
  }
  void run() {
    intra->start();
    self->start();
    self->wait();
    intra->wait();
  }
};


template <typename T>
class Gather_striped {

  int myid;
  int numproc;

  Comm::Comm<T> *self;
  Comm::Comm<T> *intra;
  Comm::Comm<T> *translation;
  Comm::Comm<T> *assembly;

  T *tempbuf;

  public:

  Gather_striped(T *sendbuf, T *recvbuf, size_t count, int root, int groupsize, const MPI_Comm &comm_mpi) {
    MPI_Comm_rank(comm_mpi, &myid);
    MPI_Comm_size(comm_mpi, &numproc);

    self = new Comm::Comm<T>(comm_mpi, Comm::IPC);
    intra = new Comm::Comm<T>(comm_mpi, Comm::IPC);
    translation = new Comm::Comm<T>(comm_mpi, Comm::NCCL);
    assembly = new Comm::Comm<T>(comm_mpi, Comm::IPC);

    int numgroup = numproc / groupsize;
    int mygroup = myid / groupsize;
    int rootgroup = root / groupsize;

    // SELF
    if(myid == ROOT) printf("self ");
    self->add(sendbuf, 0, recvbuf, root * count, count, root, root);
    //INTRA
    for(int group = 0; group < numgroup; group++)
      for(int p = 0; p < groupsize; p++) {
        int sender = group * groupsize + p;
        if(group == rootgroup) {
          if(sender != root) {
            if(myid == ROOT) printf("intra ");
            intra->add(sendbuf, 0, recvbuf, sender * count, count, sender, root);
          }
        }
      }
    // TRANSLATE
    if(mygroup == rootgroup)
      if(myid != root)
        cudaMalloc(&tempbuf, count * sizeof(T) * (numgroup - 1));
    for(int p = 0; p < groupsize; p++) {
      int tempcount = 0;
      for(int group = 0; group < numgroup; group++)
        if(group != rootgroup) {
          int sender = group * groupsize + p;
          int recver = rootgroup * groupsize + p;
          if(recver != root) {
            if(myid == ROOT) printf("translation buffer ");
            translation->add(sendbuf, 0, tempbuf, tempcount * count, count, sender, recver);
	    tempcount++;
          }
	  else {
            if(myid == ROOT) printf("translation direct ");
            translation->add(sendbuf, 0, recvbuf, sender * count, count, sender, recver);
	  }
        }
    }
    // ASSEMBLY
    for(int p = 0; p < groupsize; p++) {
      int tempcount = 0;
      for(int group = 0; group < numgroup; group++)
        if(group != rootgroup) { 
          int sender = rootgroup * groupsize + p;
          int source = group * groupsize + p;
          if(sender != root) {
            if(myid == ROOT) printf("assembly ");
            assembly->add(tempbuf, tempcount * count, recvbuf, source * count, count, sender, root);
            tempcount++;
	  }
        }
    }     
  }
  void measure_breakdown(int warmup, int numiter) {
    self->measure(warmup, numiter);
    intra->measure(warmup, numiter);
    translation->measure(warmup, numiter);
    assembly->measure(warmup, numiter);
  }
  void run() {
    translation->start();
    intra->start();
    self->start();
    self->wait();
    intra->wait();
    translation->wait();
    assembly->start();
    assembly->wait();
  }
};

template <typename T>
class Scatter_striped {

  int myid;
  int numproc;

  Comm::Comm<T> *self;
  Comm::Comm<T> *intra;
  Comm::Comm<T> *partition;
  Comm::Comm<T> *translation;

  T *tempbuf;

  public:

  Scatter_striped(T *sendbuf, T *recvbuf, size_t count, int root, int groupsize, const MPI_Comm &comm_mpi) {
    MPI_Comm_rank(comm_mpi, &myid);
    MPI_Comm_size(comm_mpi, &numproc);

    self = new Comm::Comm<T>(comm_mpi, Comm::IPC);
    intra = new Comm::Comm<T>(comm_mpi, Comm::IPC);
    partition = new Comm::Comm<T>(comm_mpi, Comm::IPC);
    translation = new Comm::Comm<T>(comm_mpi, Comm::NCCL);

    int numgroup = numproc / groupsize;
    int mygroup = myid / groupsize;
    int rootgroup = root / groupsize;

    // SELF
    if(myid == ROOT) printf("self ");
    self->add(sendbuf, root * count, recvbuf, 0, count, root, root);
    //INTRA
    for(int group = 0; group < numgroup; group++)
      for(int p = 0; p < groupsize; p++) {
        int recver = group * groupsize + p;
        if(group == rootgroup) {
          if(recver != root) {
            if(myid == ROOT) printf("intra ");
            intra->add(sendbuf, recver * count, recvbuf, 0, count, root, recver);
          }
        }
      }
    // PARTITION
    if(mygroup == rootgroup)
      if(myid != root)
        cudaMalloc(&tempbuf, count * sizeof(T) * (numgroup - 1));
    for(int p = 0; p < groupsize; p++) {
      int tempcount = 0;
      for(int group = 0; group < numgroup; group++)
        if(group != rootgroup) {
          int recver = rootgroup * groupsize + p;
          int dest = group * groupsize + p;
          if(recver != root) {
            if(myid == ROOT) printf("partition ");
            partition->add(sendbuf, dest * count, tempbuf, tempcount * count, count, root, recver);
            tempcount++;
          }
        }
    }
    // TRANSLATE
    for(int p = 0; p < groupsize; p++) {
      int tempcount = 0;
      for(int group = 0; group < numgroup; group++)
        if(group != rootgroup) {
          int recver = group * groupsize + p;
          int sender = rootgroup * groupsize + p;
          if(sender != root) {
            if(myid == ROOT) printf("translation buffer ");
            translation->add(tempbuf, tempcount * count, recvbuf, 0, count, sender, recver);
            tempcount++;
          }
          else {
            if(myid == ROOT) printf("translation direct ");
            translation->add(sendbuf, recver * count, recvbuf, 0, count, sender, recver);
          }
        }
    }
  }
  void measure_breakdown(int warmup, int numiter) {
    self->measure(warmup, numiter);
    intra->measure(warmup, numiter);
    partition->measure(warmup, numiter);
    translation->measure(warmup, numiter);
  }
  void run() {
    partition->start();
    partition->wait();
    translation->start();
    intra->start();
    self->start();
    self->wait();
    intra->wait();
    translation->wait();
  }
};

template <typename T, class Coll>
void measure(size_t count, int warmup, int numiter, Coll &coll) {

  int myid;
  int numproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);

  double times[numiter];
  if(myid == ROOT)
    printf("%d warmup iterations (in order):\n", warmup);
  for (int iter = -warmup; iter < numiter; iter++) {
    MPI_Barrier(MPI_COMM_WORLD);
    double time = MPI_Wtime();
    coll.run();
    time = MPI_Wtime() - time;
    MPI_Allreduce(MPI_IN_PLACE, &time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if(iter < 0) {
      if(myid == ROOT)
        printf("warmup: %e\n", time);
    }
    else
      times[iter] = time;
  }
  std::sort(times, times + numiter,  [](const double & a, const double & b) -> bool {return a < b;});

  if(myid == ROOT) {
    printf("%d measurement iterations (sorted):\n", numiter);
    for(int iter = 0; iter < numiter; iter++) {
      printf("time: %.4e", times[iter]);
      if(iter == 0)
        printf(" -> min\n");
      else if(iter == numiter / 2)
        printf(" -> median\n");
      else if(iter == numiter - 1)
        printf(" -> max\n");
      else
        printf("\n");
    }
    printf("\n");
    double minTime = times[0];
    double medTime = times[numiter / 2];
    double maxTime = times[numiter - 1];
    double avgTime = 0;
    for(int iter = 0; iter < numiter; iter++)
      avgTime += times[iter];
    avgTime /= numiter;
    double data = count * sizeof(T) * numproc;
          if (data < 1e3)
        printf("data: %d bytes\n", (int)data);
      else if (data < 1e6)
        printf("data: %.4f KB\n", data / 1e3);
      else if (data < 1e9)
        printf("data: %.4f MB\n", data / 1e6);
      else if (data < 1e12)
        printf("data: %.4f GB\n", data / 1e9);
      else
        printf("data: %.4f TB\n", data / 1e12);
    printf("minTime: %.4e us, %.4e s/GB, %.4e GB/s\n", minTime * 1e6, minTime / data * 1e9, data / minTime / 1e9);
    printf("medTime: %.4e us, %.4e s/GB, %.4e GB/s\n", medTime * 1e6, medTime / data * 1e9, data / medTime / 1e9);
    printf("maxTime: %.4e us, %.4e s/GB, %.4e GB/s\n", maxTime * 1e6, maxTime / data * 1e9, data / maxTime / 1e9);
    printf("avgTime: %.4e us, %.4e s/GB, %.4e GB/s\n", avgTime * 1e6, avgTime / data * 1e9, data / avgTime / 1e9);
    printf("\n");
  }
}


template <class Coll>
void validate(int *sendbuf_d, int *recvbuf_d, size_t count, int pattern, Coll &coll) {

  int myid;
  int numproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);

  int *recvbuf = new int[count * numproc];
  int *sendbuf = new int[count * numproc];

  switch(pattern) {
    case 1:
      if(myid == ROOT)
        printf("VERIFY GATHER\n");
      for(size_t i = 0; i < count; i++)
        sendbuf[i] = myid;
      cudaMemcpy(sendbuf_d, sendbuf, count * sizeof(int), cudaMemcpyHostToDevice);
      if(myid == ROOT)
        memset(recvbuf, -1, count * sizeof(int) * numproc);

      coll.run();

      if(myid == ROOT) {
        cudaMemcpy(recvbuf, recvbuf_d, count * sizeof(int) * numproc, cudaMemcpyDeviceToHost);
        bool pass = true;
        for(int p = 0; p < numproc; p++)
          for(size_t i = p * count; i < (p + 1) * count; i++) {
            if(recvbuf[i] != p)
              pass = false;
        }
        if(pass)
          printf("PASSED!\n");
        else
          printf("FAILED!!!\n");
        printf("\n");
      }
      break;
    case 2:
      if(myid == ROOT)
        printf("VERIFY SCATTER\n");
      if(myid == ROOT) {
        for(int p = 0; p < numproc; p++)
          for(size_t i = p * count; i < (p + 1) * count; i++)
            sendbuf[i] = p;
        cudaMemcpy(sendbuf_d, sendbuf, count * sizeof(int) * numproc, cudaMemcpyHostToDevice);
      }
      memset(recvbuf, -1, count * sizeof(int));


      coll.run();

     // cudaStream_t stream;
     // cudaStreamCreate(&stream);
     // cudaMemcpyAsync(recvbuf, recvbuf_d, count * sizeof(int), cudaMemcpyDeviceToHost, stream);
     //  cudaStreamSynchronize(stream);
      bool pass = true;
      for(size_t i = 0; i < count; i++) {
        // printf("myid %d recvbuf[%d] = %d\n", myid, i, recvbuf[i]);
        if(recvbuf[i] != myid)
          pass = false;
      }
      MPI_Allreduce(MPI_IN_PLACE, &pass, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);
      if(myid == ROOT) {
        if(pass)
          printf("PASSED!\n");
        else
          printf("FAILED!!!\n");
        printf("\n");
      }
    break;
  }

  delete[] sendbuf;
  delete[] recvbuf;
};

