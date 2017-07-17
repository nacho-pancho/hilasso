#ifndef MDLS_HUBER
#define MDLS_HUBER

template<typename T>
struct sc_huber_workspace {

  sc_huber_workspace(): Hk(),Rj(),Gk(),Yk(),YtHk() {}
  
  sc_huber_workspace(int M) {
    this->allocate(M);
  }

  void allocate(int M) {
    Hk.allocate(M,M);
    Gk.allocate(M);
    Rj.allocate(M);
    Yk.allocate(M);
    YtHk.allocate(M);
    prevDk.allocate(M);
  }
  
  void free() {
    if (Hk.data != NULL) {
      Hk.free();
    }
    if (Gk.data != NULL)  {
      Gk.free();
    }
    if (Rj.data != NULL) {
      Rj.free();
    }
    if (Yk.data != NULL) {
      Yk.free();
    }
    if (YtHk.data != NULL) {
      YtHk.free();
    }
  }

  ~sc_huber_workspace() {
    this->free();
  }
  //
  //
  //
  sc_matrix<T> Hk; // inverse of Hessian matrix for atom k
  sc_vector<T> Rj;    // residual for sample j
  sc_vector<T> Gk;    // gradient for atom k
  //
  // BFGS auxiliary variables
  //
  sc_vector<T> Yk;   // Gk(t)-Gk(t-1) at iteration t > 1
  sc_vector<T> YtHk; // Yk^T * inv(H_k)
  sc_vector<T> prevDk;
};
#endif
