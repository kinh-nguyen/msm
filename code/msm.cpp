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

// Log-skew-logistic baseline hazard function 
template <class Type>
Type LSLx(Type t, Type lamda, Type p, Type gamma) {
  if (t == 0.0) return Type(0.0); // can be undefined without reparameterization
  Type
      u = pow(lamda * t, -p),
      n = gamma * p * u,
      d = t * (1 + u) * (pow(1 + u, gamma) - Type(1.0));
  return n / d;
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

  // log-skew-logistic hazard model
  PARAMETER_VECTOR(log_lambda)
  vector<Type> lambda = exp(log_lambda);
  dll -= dnorm(lambda, Type(0), Type(1), true).sum() + log_lambda.sum();

  PARAMETER_VECTOR(log_p);
  vector<Type> p = exp(log_p);
  dll -= dnorm(p, Type(0), Type(1), true).sum() + log_p.sum();

  PARAMETER_VECTOR(log_gamma);
  vector<Type> gamma = exp(log_gamma);
  dll -= dnorm(gamma, Type(0), Type(1), true).sum() + log_gamma.sum();

  PQ<Type> KM;
  vector<Type> q_rs(N_PAR);

  // Add vector to track individual log-likelihoods
  vector<Type> indiv_ll(A.size());
  indiv_ll.setZero();

  for (int i = 0; i < A.size(); i++) {
    Type ll = 0;
    for (int j = 0; j < N_PAR; j++)
      q_rs[j] = LSLx(Type(start[i]), lambda[j], p[j], gamma[j]); 
    if (fit[i] == 0) 
      ll = -n[i] * KM(q_rs, true, true)(A[i], Z[i]); 
    else if (fit[i] == 1) 
      ll = -n[i] * KM(q_rs, true, false)(A[i], Z[i]); 
    else if (fit[i] == 2) {
      matrix<Type> cumP = matrix<Type>::Identity(N_Q, N_Q);
      for (int t = start[i]; t < end[i]; t++)
      {
        for (int j = 0; j < N_PAR; j++) 
          q_rs[j] = LSLx(Type(t), lambda[j], p[j], gamma[j]);
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
        eta_rep[j] = LSLx(Type(s), lambda[j], p[j], gamma[j]);
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
    REPORT(p);
    REPORT(gamma);
    REPORT(lambda);
    // Report individual log-likelihoods
    REPORT(indiv_ll);
  }
  return dll;
}
