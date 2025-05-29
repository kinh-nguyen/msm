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
  parallel_accumulator<Type> dll(this);

  // data
  DATA_IVECTOR(A);
  DATA_IVECTOR(Z);
  DATA_IVECTOR(start);
  DATA_IVECTOR(end);
  DATA_VECTOR(n);
  DATA_IVECTOR(fit);
  DATA_INTEGER(n_age);
  // for padding non-data 
  DATA_VECTOR(age_dv);
  DATA_IVECTOR(minage);
  DATA_IVECTOR(n_epis);
  // priors
  DATA_VECTOR(prior_base);

  PARAMETER_VECTOR(intercepts);
  dll -= dnorm(intercepts, prior_base(0), prior_base(1), true).sum();

  // Age AR(1) model
  PARAMETER_VECTOR(pacf_vec); // length 1 * N_PAR
  PARAMETER_VECTOR(age_sm); // length depends on minage and n_epis
  dll -= dnorm(pacf_vec, Type(0), Type(1), true).sum();
  int pid = 0, did = 0;
  for (int i = 0; i < N_PAR; i++) {
    Type phi = 2. * exp(pacf_vec[i]) / (1. + exp(pacf_vec[i])) - 1.;
    vector<Type> age_sm_ = age_sm(seqN(pid, n_epis[i]));
    dll += density::AR1(phi)(age_sm_);
    age_dv(seqN(did + minage[i], n_epis[i])) = age_sm_; // "padding" with zeros
    did += n_age; // advance to next transition in padded vector
    pid += n_epis[i]; // advance to next parameter in estimated parameter vector
  }

  PQ<Type> KM;
  vector<Type> eta(N_PAR);

  for (int i = 0; i < A.size(); i++) {
    for (int j = 0; j < N_PAR; j++) 
      eta[j] = exp(intercepts[j] + age_dv(j * n_age + start[i]));
    if (fit[i] == 0) 
      dll -= n[i] * KM(eta, true, true)(A[i], Z[i]); 
    else if (fit[i] == 1) 
      dll -= n[i] * KM(eta, true, false)(A[i], Z[i]); 
    if (fit[i] == 2) {
      matrix<Type> cumP = matrix<Type>::Identity(N_Q, N_Q);
      for (int t = start[i]; t < end[i]; t++)
      {
        for (int j = 0; j < N_PAR; j++) 
          eta[j] = exp(intercepts[j] + age_dv(j * n_age + t));
        cumP = cumP * KM(eta, false, false);
      }
      dll -= n[i] * log(cumP(A[i], Z[i]));
    } 
  }
  
  SIMULATE {
    array<Type> qM_rep(N_Q, N_Q, n_age), pM_rep(N_Q, N_Q, n_age);
    for (int s = 0; s < n_age; ++s) {
      vector<Type> eta_rep(N_PAR);
      for (int j = 0; j < N_PAR; j++)
        eta_rep[j] = exp(intercepts[j] + age_dv(j * n_age + s));
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
    REPORT(intercepts);
    REPORT(pacf_vec);
    REPORT(age_dv);
    REPORT(age_sm);
  }
  return dll;
}
