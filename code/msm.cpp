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
      if (isQ)
        return isLog ? qM.array().log().matrix() : qM;
      pM = expm(qM);
      return isLog ? pM.array().log().matrix() : pM;
    };
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  Type dll = 0.0;
  Type prior = 0.0;

  // data
  DATA_IVECTOR(A);
  DATA_IVECTOR(Z);
  DATA_VECTOR(start);
  DATA_VECTOR(end);
  DATA_VECTOR(n);
  DATA_IVECTOR(cens); // censored indicator
  
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
    for (int j = 0; j < N_PAR; j++) 
      eta[j] = exp(intercepts[j] + betas[j] * start[i]); 
    if (cens[i] == 0 && (A[i] - Z[i] != 0)) 
      ll[i] = n[i] * KM(eta, start[i], end[i], true, true)(A[i], Z[i]); 
    else 
      ll[i] = n[i] * KM(eta, start[i], end[i], true, false)(A[i], Z[i]); 
    dll -= ll[i];
  }
  dll += prior;
  REPORT(intercepts);
  REPORT(betas);
  REPORT(ll);
  array<Type> qM_rep(N_Q, N_Q, 50), pM_rep(N_Q, N_Q, 50);
  for (int s = 0; s <= 49; ++s) {
    vector<Type> eta_rep(N_PAR);
    for (int j = 0; j < N_PAR; j++)
      eta_rep[j] = exp(intercepts[j] + betas[j] * s);
    qM_rep.col(s) = Eigen::Map<Matrix<Type, N_Q, N_Q> >(KM(eta_rep, s, s + 1, false, true).data()).reshaped(N_Q * N_Q, 1);
    pM_rep.col(s) = Eigen::Map<Matrix<Type, N_Q, N_Q> >(KM(eta_rep, s, s + 1, false, false).data()).reshaped(N_Q * N_Q, 1);
  }
  REPORT(qM_rep);
  REPORT(pM_rep);
  REPORT(qM_rep);
  REPORT(pM_rep);
  return dll;
}
