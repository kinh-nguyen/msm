#include <TMB.hpp>
#include "ktools.hpp"

#define N_PAR 7
#define A_REF 20 // referenced age 20
#define N_Q 7

using Eigen::seqN;
using ktools::rw1_nll;

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
  DATA_IVECTOR(len_dv);
  DATA_IVECTOR(len_lv);
  DATA_IVECTOR(len_pv);
  DATA_IVECTOR(minage);
  DATA_IVECTOR(maxage);
  
  // priors
  DATA_VECTOR(prior_base);

  PARAMETER_VECTOR(intercepts);
  dll -= dnorm(intercepts, prior_base(0), prior_base(1), true).sum();

  // Age RW(1) model
  PARAMETER(log_sigma_rw1); // share variance across transitions
  PARAMETER_VECTOR(age_sm); // length = sum(len_pv)

  vector<Type> age_dv(len_dv.sum());
  age_dv.setZero(); 

  int pid = 0, did = 0;
  for (int i = 0; i < N_PAR; i++) {
    vector<Type> age_pv = age_sm(seqN(pid, len_pv[i]));
    vector<Type> age_lv(len_lv[i]);
    age_lv.setZero(); 
    int s1 = A_REF - minage[i]; 
    int s2 = maxage[i] - A_REF; 
    age_lv(seqN(0, s1)) = age_pv(seqN(0, s1)); // left of A_REF
    age_lv(seqN(s1 + 1, s2)) = age_pv(seqN(s1, s2)); // right of A_REF
    age_dv(seqN(did + minage[i], len_lv[i])) = age_lv; // store for linear predictor
    dll += rw1_nll(age_lv, log_sigma_rw1, Type(0.0), Type(0.5), true, false);
    pid += len_pv[i]; 
    did += len_dv[i]; 
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
    REPORT(log_sigma_rw1);
    REPORT(age_dv);
    REPORT(age_sm);
  }
  return dll;
}
