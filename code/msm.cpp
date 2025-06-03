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

// https://dlmf.nist.gov/8.2
// t is time, exp_b0 = alpha, b1 is beta > 0
template <class Type>
Type IGx(Type t, Type log_t, Type a, Type b) {
  Type log_igm = lgamma(a) + pgamma(b / t, a, Type(1.0));
  Type lhz = a * log(b) - (a + Type(1.0)) * log_t - b / t - log_igm;
  return exp(lhz);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);

  // data
  DATA_IVECTOR(A);
  DATA_IVECTOR(Z);
  DATA_IVECTOR(istart);
  DATA_VECTOR(dstart);
  vector<Type> log_t = log(dstart);
  DATA_IVECTOR(end);
  DATA_VECTOR(n);
  DATA_IVECTOR(fit);
  DATA_INTEGER(n_age);

  PARAMETER_VECTOR(log_mu);
  DATA_VECTOR(prior_mu);
  vector<Type> mu = exp(log_mu);
  dll -= dnorm(log_mu, prior_mu(0), prior_mu(1), true).sum();

  PARAMETER_VECTOR(log_sigma);
  DATA_VECTOR(prior_sigma);
  vector<Type> sigma = exp(log_sigma);
  dll -= dnorm(log_sigma, prior_sigma(0), prior_sigma(1), true).sum();

  // Transformation
  vector<Type> a = Type(1.0) / (sigma * sigma);
  vector<Type> b = mu * (a + Type(1.0));

  PQ<Type> KM;
  vector<Type> q_rs(N_PAR);

  vector<Type> indiv_ll(A.size());
  indiv_ll.setZero();

  auto fill_qrs = [&](Type age, Type log_age) {
    for (int p = 0; p < N_PAR; p++)
      q_rs[p] = IGx(age, log_age, a[p], b[p]);
  };

  for (int i = 0; i < A.size(); i++) {
    Type ll = 0;
    fill_qrs(dstart[i], log_t[i]);
    if (fit[i] == 0) 
      ll = -n[i] * KM(q_rs, true, true)(A[i], Z[i]); 
    else if (fit[i] == 1) 
      ll = -n[i] * KM(q_rs, true, false)(A[i], Z[i]); 
    else if (fit[i] == 2) {
      matrix<Type> cumP = matrix<Type>::Identity(N_Q, N_Q);
      for (int t = istart[i]; t < end[i]; t++)
      {
        fill_qrs(Type(t), log(Type(t)));
        cumP = cumP * KM(q_rs, false, false);
      }
      ll = -n[i] * log(cumP(A[i], Z[i]));
    }
    dll += ll;
    indiv_ll[i] = ll;
  }
  
  SIMULATE {
    array<Type> qM_rep(N_Q, N_Q, n_age), pM_rep(N_Q, N_Q, n_age);
    for (int s = 1; s < n_age; ++s) {
      fill_qrs(Type(s), log(Type(s)));
      matrix<Type> qM_mat = KM(q_rs, false, true);
      matrix<Type> pM_mat = KM(q_rs, false, false);
      for (int i = 0; i < N_Q; ++i)
        for (int j = 0; j < N_Q; ++j) {
          qM_rep(i, j, s) = qM_mat(i, j);
          pM_rep(i, j, s) = pM_mat(i, j);
        }
    }
    REPORT(qM_rep);
    REPORT(pM_rep);
    REPORT(mu);
    REPORT(sigma);
    REPORT(indiv_ll);
  }
  return dll;
}
