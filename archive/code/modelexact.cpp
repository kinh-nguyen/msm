#include <TMB.hpp>
#include "ktools.hpp"

#define N_PAR 7
#define N_Q 7

using Eigen::seqN;

template <class T> 
struct PQ {
    matrix<T> 
        qM = matrix<T>(N_Q, N_Q),
        pM = matrix<T>(N_Q, N_Q);
    PQ () {};
    matrix<T> operator() (vector<T> eta, T start, T end, bool isLog = true, bool isQ = false) {
      qM.setZero(); pM.setZero();
      qM(0, 1) = eta(0); // debut
      qM(0, 2) = eta(1); // marriage from virgin
      qM(1, 2) = eta(2); // marriage from debut
      qM(2, {3,4,5}) = eta({3,4,5}); // marriage dissolution
      qM({3,4,5}, 6) = eta({6,6,6}); // disso > remarried
      qM(6, {3,4,5}) = eta({3,4,5}); // remarried > disso = married > disso, we could add a(three) scaling parameter as well?
      qM.diagonal() = T(-1) * qM.rowwise().sum();
      qM *= (end - start);
      pM = expm(qM);
      if (isQ) {
        if (isLog) qM = log(qM.array());
        return qM;
      } 
      if (isLog) pM = log(pM.array());
      return pM;
    };
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  Type dll = 0.0;
  // parallel_accumulator<Type> dll(this);
  Type prior = 0.0;

  // data
  DATA_IVECTOR(A);
  DATA_IVECTOR(Z);
  DATA_VECTOR(start);
  DATA_VECTOR(end);
  DATA_VECTOR(n);
  
  vector<Type> age_start = start/100;
  vector<Type> age_end = end/100;

  // priors
  DATA_VECTOR(prior_base);

  PARAMETER_VECTOR(intercepts);
  prior -= dnorm(intercepts, prior_base(0), prior_base(1), true).sum();

  PARAMETER_VECTOR(betas);
  prior -= dnorm(betas, prior_base(0), prior_base(1), true).sum();

  PQ<Type> KM;
  vector<Type> eta(N_PAR);
  vector<Type> ll(A.size());

  for (int i = 0; i < A.size(); i++) {
    if (A[i] == Z[i]) { // Exact, remain
      for (int j = 0; j < N_PAR; j++) 
        eta[j] = exp(intercepts[j] + betas[j] * age_start[i]); // 
      ll[i] = n[i] * KM(eta, start[i], end[i], false, true)(A[i], Z[i]); // take Q only
    }
    else if (A[i] == 2 && Z[i] > 2) { // Censoring transition, known end state
      for (int j = 0; j < N_PAR; j++) 
        eta[j] = exp(intercepts[j] + betas[j] * age_start[i]); //aam
      ll[i] = n[i] * KM(eta, start[i], end[i], true, false)(A[i], Z[i]); // only one with P instead of T
    } 
    else { // Exact, move
      for (int j = 0; j < N_PAR; j++) 
        eta[j] = exp(intercepts[j] + betas[j] * age_start[i]); // age
      ll[i] += n[i] * KM(eta, Type(0), Type(1), true, true)(A[i],Z[i]); // take log Q01 x 1
    }
    dll -= ll[i];
  }
  dll += prior;
  REPORT(intercepts);
  REPORT(betas);
  REPORT(ll);
  return dll;
}
