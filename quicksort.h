#include <algorithm>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <random>

#include "parallel.h"


using namespace parlay;

//hyper params, set to 1e5 and false before submitting
size_t threshhold = 300000; 
size_t threshhold2 = 1000000;
bool print_ = false;

template <class T>
void print(T* A, size_t n){
  for(size_t i = 0; i < n; ++i){
    std::cout << A[i] << " ";
  }
  std::cout<<std::endl;
}

inline uint64_t hash64_q(uint64_t u) {
  uint64_t v = u * 3935559000370003845ul + 2691343689449507681ul;
  v ^= v >> 21;
  v ^= v << 37;
  v ^= v >> 4;
  v *= 4768777513237032717ul;
  v ^= v << 20;
  v ^= v >> 41;
  v ^= v << 5;
  return v;
}

template <class T>
size_t Sum(T* A,size_t n){
  //Right now it is linear, change to parlallel scan
  for(size_t i = 1; i < n; ++i){
    A[i] += A[i-1];
  }
  return A[n-1];
}

template<class T>
T reduce(T* A, size_t n){
  if (n < threshhold2){
    T sum = 0;
    for(size_t i = 0; i < n; ++i){
      sum+=A[i];
    }
    return sum;
  }
  else{
    T v1,v2;
    auto f1 = [&]() {v1 = reduce(A,n/2);};
    auto f2 = [&]() {v2 = reduce(A+n/2,n-n/2);};
    par_do(f1,f2);
    return v1 + v2;
  }
}

template <class T>
void scan_up(T* A, size_t s, size_t t){
  if (s == t){
    return;
  }
  else{
    size_t m = (s+t)/2;
    if (s+t < threshhold2){
      scan_up(A,s,m);
      scan_up(A,m+1,t);
    }
    else{
      auto f1 = [&]() {scan_up(A,s,m);};
      auto f2 = [&]() {scan_up(A,m+1,t);};
      par_do(f1,f2);
    }
    A[t] += A[m];
  }
}

template <class T>
void scan_down(T* A,size_t s, size_t t, T p){
  if(print_){std::cout<< "(s,t,p,ls): (" << s << "," << t << "," << p << "," << A[(s+t)/2] <<  ")" << std::endl;}
  if (s == t){
    A[s] = p;
  }
  else{
    size_t m = (s+t)/2;
    T leftsum = A[m];
    if (s+t < threshhold){
      scan_down(A,s,m,p);
      scan_down(A,m+1,t,leftsum+p);
    }
    else{
      auto f1 = [&]() {scan_down(A,s,m,p);};
      auto f2 = [&]() {scan_down(A,m+1,t,leftsum+p);};
      par_do(f1,f2);
    }
  }
  return;
}

template <class T>
size_t Sum_parallel(T* A, T* B,size_t n){
  size_t k = std::sqrt(n);
  size_t j = n/k;
  if(print_){std::cout << "(n,k): " << '(' << n << ',' << k << ')' << std::endl;}
  //Compute sums for each element
  B[j-1] = reduce(A,j);
  parallel_for(2,k+1,[&](size_t i){
    B[std::min(n,i*j)-1] = reduce(A+((i-1)*j),j);
  });

  //Sums computed, add sums linearly
  for(size_t i = 2; i <= k; ++i){
    B[std::min(n,i*j)-1] += B[(i-1)*j-1];
  }

  //Compute prefix sum
  parallel_for(0,k+1,[&](size_t i){
    for(size_t l = 0; l < j; ++l){
      if (std::min(i*j + l,n-1) == 0){
        B[std::min(i*j + l,n-1)] = A[std::min(i*j + l,n-1)];
      }
      else{
        B[std::min(i*j + l,n-1)] = B[std::min(i*j + l,n-1)-1]+A[std::min(i*j + l,n-1)];
      }
    }
  });
  if(print_){std::cout<<B[n-1]<<std::endl;}
  return B[n-1];
}

template <class T>
void filter(T *A, int* B, T* C, size_t n){
  if (B[0] == 1){
    C[0] = A[0];
  }
  parallel_for(1,n,[&](size_t i){
    if (B[i]!=B[i-1]){
      C[B[i]-1] = A[i];
    }
  });
}

template<class T>
void filter_back(T*A,int*B,T*C,size_t n){
  if (B[0] == 1){
    C[n-1] = A[0];
  }
  parallel_for(1,n,[&](size_t i){
    if (B[i] != B[i-1]){
      C[n-B[i]] = A[i];
    }
  });
}

template<class T>
void flags_less(int *B, T*A, T* C, size_t n, T pivot){
  parallel_for(0,n,[&](size_t i){
   });
}

template<class T>
void flags_eq(int *B, T*A, size_t n, T pivot){
  parallel_for(0,n,[&](size_t i){
    B[i] = (A[i] == pivot);
  });
}

template<class T>
void flags_great(int *B, T*A, size_t n, T pivot){
  parallel_for(0,n,[&](size_t i){
    B[i] = (A[i] > pivot);
  });
}

template <class T>
void flags_both(int*B, T*A, T* C, int* E, size_t n, T pivot){
  parallel_for(0,n,[&](size_t i){
    B[i] = (A[i] < pivot);
    E[i] = (A[i] > pivot);
    C[i] = A[i];
  });
}

template <class T>
void filter_both(T *A, int* B, T* C, int* D, size_t n){ //A has elements, put them into C... B has lesser, D has greater
  if (B[0] == 1){
    C[0] = A[0];
  }
  if(D[0] == 1){
    C[n-1] = A[0];
  }
  parallel_for(1,n,[&](size_t i){
    if (B[i]!=B[i-1]){
      C[B[i]-1] = A[i];
    }
    if (D[i] != D[i-1]){
      C[n-D[i]] = A[i];
    }
  });
}

template <class T>
void partition(T *A, int* B, T* C, int* D, int* E, size_t n, T pivot,size_t* mid1, size_t* mid2){
  flags_both(B,A,C,E,n,pivot); //C contains A, B contains flags for less than, E contains flags for greater than. Trash: A
  *mid1 = Sum_parallel(B,D,n); //D contains prefix sum of B. Trash: A,B
  *mid2 = Sum_parallel(E,B,n); //B contains prefix sum of E. Trash: E,A
  filter_both(C,D,A,B,n); //C has been put back into A, 
  parallel_for(*mid1,n-(*mid2),[&](size_t i){ //Adding pivot occurances
    A[i] = pivot;
  });

  *mid2 = (n - *mid2);
  return;
}

template <class T>
void quicksort_r(T *A, int* B, T* C, int* D, int* E, size_t n){
  if (n <= threshhold){
    std::sort(A,A+n);
  }
  else{
    T pivots [1000];
    for(int i = 0 ;i < 1000 ; i++){
      pivots[i] = A[hash64_q(i*12345)%n];
    }
    std::sort(std::begin(pivots),std::end(pivots));

    //std::pair<size_t,size_t> mids;
    size_t mid1;
    size_t mid2;
    partition(A,B,C,D,E,n,pivots[500],&mid1,&mid2);
    //size_t mid1 = mids.first;
    //size_t mid2 = mids.second;
    auto f1 = [&]() {quicksort_r(A,B,C,D,E,mid1);};
    auto f2 = [&]() {quicksort_r(A+mid2,B+mid2,C+mid2,D+mid2,E+mid2,n-mid2);};
    par_do(f1,f2);
    if(print_){std::cout << "Mid: " << mid1 << " " << mid2 << std::endl;}
  }
  return;
}


template <class T>
void quicksort(T *A, size_t n) {
  /*
  std::mt19937 mt_rng();
  std::cout << mt_rng;
  return;*/
  if(print_){std::cout << "Initial array: "; print(A,n);}
  int* B = (int*)malloc(n * sizeof(int));
  T* C = (T*)malloc(n * sizeof(T));
  int* D = (int*)malloc(n * sizeof(int));
  int* E = (int*)malloc(n * sizeof(int));
  quicksort_r(A,B,C,D,E,n);
  if(print_){std::cout << "Initial array: "; print(A,n);}
  if(print_){std::cout << "Initial array: "; print(C,n);}
  return;
}
