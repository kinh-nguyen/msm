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
    matrix<T> operator() (vector<T> q_rs, bool isLog = true, bool isQ = false) {
      qM.setZero(); pM.setZero();
      qM(0, 1) = q_rs(0); // debut
      qM(0, 2) = q_rs(1); // marriage from virgin
      qM(1, 2) = q_rs(2); // marriage from debut
      qM(2, {3,4,5}) = q_rs({3,4,5}); // marriage dissolution
      qM({3,4,5}, 6) = q_rs({6,6,6}); // disso > remarried
      qM(6, {3,4,5}) = q_rs({3,4,5}); // remarried > disso = married > disso, we could add a(three) scaling parameter as well?
      qM.diagonal() = T(-1) * qM.rowwise().sum();
      if (isQ)
        return isLog ? qM.array().log().matrix() : qM;
      pM = atomic::expm(qM);
      return isLog ? pM.array().log().matrix() : pM;
    };
};

// Log-logistic baseline hazard function where
// t is time, exp_b0 = alpha, b1 is beta > 0
template <class Type>
Type LLGx(Type t, Type exp_b0, Type b1) {
  if (t == 0.0) return Type(0.0); 
  Type 
    t_pow_b1 = pow(t, b1),
    numerator = exp_b0 * b1 * t_pow_b1 / t,
    denominator = Type(1.0) + exp_b0 * t_pow_b1;
  return numerator / denominator;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);

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

  // log-logistic hazard model
  PARAMETER_VECTOR(b0);
  dll -= dnorm(b0, prior_base(0), prior_base(1), true).sum();
  vector<Type> exp_b0 = exp(b0);

  PARAMETER_VECTOR(log_b1)
  vector<Type> b1 = exp(log_b1);
  dll -= dnorm(b1, Type(0), Type(1), true).sum() + log_b1.sum(); 

  PQ<Type> KM;
  vector<Type> q_rs(N_PAR);

  // Add vector to track individual log-likelihoods
  vector<Type> indiv_ll(A.size());
  indiv_ll.setZero();

  for (int i = 0; i < A.size(); i++) {
    Type ll = 0;
    for (int j = 0; j < N_PAR; j++)
      q_rs[j] = LLGx(Type(start[i]), exp_b0[j], b1[j]); 
    if (fit[i] == 0) 
      ll = -n[i] * KM(q_rs, true, true)(A[i], Z[i]); 
    else if (fit[i] == 1) 
      ll = -n[i] * KM(q_rs, true, false)(A[i], Z[i]); 
    else if (fit[i] == 2) {
      matrix<Type> cumP = matrix<Type>::Identity(N_Q, N_Q);
      for (int t = start[i]; t < end[i]; t++)
      {
        for (int j = 0; j < N_PAR; j++) 
          q_rs[j] = LLGx(Type(t), exp_b0[j], b1[j]);
        cumP = cumP * KM(q_rs, false, false);
      }
      ll = -n[i] * log(cumP(A[i], Z[i]));
    }
    dll += ll;
    indiv_ll[i] = ll;
  }
  
  SIMULATE {
    array<Type> qM_rep(N_Q, N_Q, n_age), pM_rep(N_Q, N_Q, n_age);
    for (int s = 0; s < n_age; ++s) {
      vector<Type> eta_rep(N_PAR);
      for (int j = 0; j < N_PAR; j++)
        eta_rep[j] = LLGx(Type(s), exp_b0[j], b1[j]);
      matrix<Type> qM_mat = KM(eta_rep, false, true);
      matrix<Type> pM_mat = KM(eta_rep, false, false);
      for (int i = 0; i < N_Q; ++i)
        for (int j = 0; j < N_Q; ++j) {
          qM_rep(i, j, s) = qM_mat(i, j);
          pM_rep(i, j, s) = pM_mat(i, j);
        }
    }
    REPORT(qM_rep);
    REPORT(pM_rep);
    REPORT(b0);
    REPORT(b1);
    // Report individual log-likelihoods
    REPORT(indiv_ll);
  }
  return dll;
}
