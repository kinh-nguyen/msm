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
    matrix<T> operator() (vector<T> eta, bool isLog = true, bool isQ = false) {
      qM.setZero(); pM.setZero();
      qM(0, 1) = eta(0); // debut
      qM(0, 2) = eta(1); // marriage from virgin
      qM(1, 2) = eta(2); // marriage from debut
      qM(2, {3,4,5}) = eta({3,4,5}); // marriage dissolution
      qM({3,4,5}, 6) = eta({6,6,6}); // disso > remarried
      qM(6, {3,4,5}) = eta({3,4,5}); // remarried > disso = married > disso, we could add a(three) scaling parameter as well?
      qM.diagonal() = T(-1) * qM.rowwise().sum();
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
  DATA_IVECTOR(start);
  DATA_IVECTOR(end);
  DATA_VECTOR(n);
  DATA_IVECTOR(fit);
  DATA_INTEGER(n_age);
  
  // priors
  DATA_VECTOR(prior_base);

  PARAMETER_VECTOR(intercepts);
  prior -= dnorm(intercepts, prior_base(0), prior_base(1), true).sum();

  // Age ARk
  PARAMETER_VECTOR(pacf_vec); // length 2 * N_PAR
  PARAMETER_VECTOR(age_sm); // length N_PAR * n_age (n_age = 50)
  for (int i = 0; i < N_PAR; i++) {
    vector<Type> age_sm_ = age_sm(seqN(i * n_age, n_age));
    vector<Type> pacf_ = pacf_vec(seqN(i * 2, 2));
    dll -= dnorm(age_sm_[0], Type(0.0), Type(0.001), true); 
    dll += ktools::AR2ll(pacf_, age_sm_);
  }

  PQ<Type> KM;
  vector<Type> eta(N_PAR);
  vector<Type> ll(A.size());

  for (int i = 0; i < A.size(); i++) {
    for (int j = 0; j < N_PAR; j++) 
      eta[j] = exp(intercepts[j] + age_sm(j * n_age + start[i]));
    if (fit[i] == 0) 
      ll[i] = n[i] * KM(eta, true, true)(A[i], Z[i]); 
    else if (fit[i] == 1) 
      ll[i] = n[i] * KM(eta, true, false)(A[i], Z[i]); 
    if (fit[i] == 2) {
      matrix<Type> cumP = matrix<Type>::Identity(N_Q, N_Q);
      for (int t = start[i]; t < end[i]; t++)
      {
        for (int j = 0; j < N_PAR; j++) 
          eta[j] = exp(intercepts[j] + age_sm(j * n_age + t));
        cumP = cumP * KM(eta, false, false);
      }
      ll[i] = n[i] * log(cumP(A[i], Z[i]));
    } 
  }
  dll -= ll.sum();
  dll += prior;
  REPORT(intercepts);
  REPORT(age_sm);
  REPORT(ll);
  // Convert list of matrices to 3D arrays for TMB reporting
  array<Type> qM_rep(N_Q, N_Q, n_age), pM_rep(N_Q, N_Q, n_age);
  for (int s = 0; s < n_age; ++s) {
    vector<Type> eta_rep(N_PAR);
    for (int j = 0; j < N_PAR; j++)
      eta_rep[j] = exp(intercepts[j] + age_sm(j * n_age + s));
    qM_rep.col(s) = Eigen::Map<Matrix<Type, N_Q, N_Q> >(KM(eta_rep, false, true).data()).reshaped(N_Q * N_Q, 1);
    pM_rep.col(s) = Eigen::Map<Matrix<Type, N_Q, N_Q> >(KM(eta_rep, false, false).data()).reshaped(N_Q * N_Q, 1);
  }
  REPORT(qM_rep);
  REPORT(pM_rep);
  return dll;
}
