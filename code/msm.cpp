#include <TMB.hpp>
#include "ktools.hpp"

#define N_PAR 7
#define N_Q 7
#define eps 10. * CppAD::numeric_limits<double>::epsilon()

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

template<class Type>
Type logSHASHz(Type t, Type mu, Type sigma, Type nu, Type tau) {
  Type z = (log(t) - mu) / sigma,
    shz = log(z + sqrt(z * z + 1)),
    e1 = exp(tau * shz),
    e2 = exp(-nu * shz),
    r = 0.5 * (e1 - e2),
    c = 0.5 * (tau * e1 + nu * e2),
    log_d = log(c) - 0.5 * r * r - 0.5 * log(2.0 * M_PI) - log(sigma)- 0.5 * log(1 + z * z) - log(t),
    log_s = log(1 - pnorm(r) + eps); 
  return exp(log_d - log_s);
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

  PARAMETER_VECTOR(mu);
  DATA_VECTOR(prior_mu);
  vector<Type> exp_mu = exp(mu);
  dll -= dnorm(mu, prior_mu(0), prior_mu(1), true).sum();

  PARAMETER_VECTOR(log_sigma);
  DATA_VECTOR(prior_sigma);
  vector<Type> sigma = exp(log_sigma);
  dll -= dnorm(sigma, prior_sigma(0), prior_sigma(1), true).sum() + log_sigma.sum();

  PARAMETER_VECTOR(log_nu);
  DATA_VECTOR(prior_nu);
  vector<Type> nu = exp(log_nu);
  dll -= dnorm(nu, prior_nu(0), prior_nu(1), true).sum() + log_nu.sum();
  
  PARAMETER_VECTOR(log_tau);
  DATA_VECTOR(prior_tau);
  vector<Type> tau = exp(log_tau);
  dll -= dnorm(tau, prior_tau(0), prior_tau(1), true).sum() + log_tau.sum();

  DATA_VECTOR(prior_coef);
  
  PARAMETER_VECTOR(itc);
  dll -= dnorm(itc, prior_coef(0), prior_coef(1), true).sum();
  PARAMETER_VECTOR(btt);
  dll -= dnorm(btt, prior_coef(0), prior_coef(1), true).sum();
  PARAMETER_VECTOR(tsq);
  dll -= dnorm(tsq, prior_coef(0), prior_coef(1), true).sum();

  // rw1
  DATA_MATRIX(Q);
  PARAMETER_VECTOR(m_vm);
  PARAMETER(log_qvm);
  Type qvm = exp(log_qvm);
  dll -= dnorm(qvm, Type(0), Type(1), true) + log_qvm;
  dll += ktools::rw(m_vm, Q, qvm, Type(1), true, true);

  DATA_MATRIX(Qvx);
  PARAMETER_VECTOR(m_vx);
  PARAMETER(log_qvx);
  Type qvx = exp(log_qvx);
  dll -= dnorm(qvx, Type(0), Type(1), true) + log_qvx;
  dll += ktools::rw(m_vx, Qvx, qvx, Type(1), true, true);
  
  DATA_MATRIX(Qxm);
  PARAMETER_VECTOR(m_xm);
  PARAMETER(log_qxm);
  Type qxm = exp(log_qxm);
  dll -= dnorm(qxm, Type(0), Type(1), true) + log_qxm;
  dll += ktools::rw(m_xm, Qxm, qxm, Type(1), true, true);

  PQ<Type> KM;
  vector<Type> q_rs(N_PAR);

  vector<Type> indiv_ll(A.size());
  indiv_ll.setZero();
  Type tmid = 20; // reference age

  auto fill_qrs = [&](Type age) {
  {
    Type tmid_age = age - tmid;
    int ida = CppAD::Integer(age);
    for (int p = 0; p < N_PAR; p++)
    {
      Type eta = itc[p] + btt[p] * tmid_age + tsq[p] * tmid_age * tmid_age;
      if (p == 1)
        q_rs[p] = exp(eta + m_vm[ida]);
      else if (p == 0)
        q_rs[p] = exp(eta + m_vx[ida]);
      else if (p == 2)
        q_rs[p] = exp(eta + m_xm[ida]);
      else
        q_rs[p] = logSHASHz(age, mu[p], sigma[p], nu[p], tau[p]) * exp(eta);
    }
  };

  for (int i = 0; i < A.size(); i++)
  {
    Type ll = 0;
    fill_qrs(dstart[i]);
    if (fit[i] == 0)
      ll = -n[i] * KM(q_rs, true, true)(A[i], Z[i]);
    else if (fit[i] == 1)
      ll = -n[i] * KM(q_rs, true, false)(A[i], Z[i]);
    else if (fit[i] == 2)
    {
      matrix<Type> cumP = matrix<Type>::Identity(N_Q, N_Q);
      for (int t = istart[i]; t < end[i]; t++)
      {
        fill_qrs(Type(t));
        cumP = cumP * KM(q_rs, false, false);
      }
      ll = -n[i] * log(cumP(A[i], Z[i]));
    }
    dll += ll;
    indiv_ll[i] = ll;
  }

  SIMULATE
  {
    array<Type> qM_rep(N_Q, N_Q, n_age), pM_rep(N_Q, N_Q, n_age);
    for (int s = 1; s < n_age; ++s)
    {
      fill_qrs(Type(s));
      matrix<Type> qM_mat = KM(q_rs, false, true);
      matrix<Type> pM_mat = KM(q_rs, false, false);
      for (int i = 0; i < N_Q; ++i)
        for (int j = 0; j < N_Q; ++j)
        {
          qM_rep(i, j, s) = qM_mat(i, j);
          pM_rep(i, j, s) = pM_mat(i, j);
        }
    }
    REPORT(qM_rep);
    REPORT(pM_rep);
    REPORT(exp_mu);
    REPORT(sigma);
    REPORT(tau);
    REPORT(nu);
    REPORT(m_vm);
    REPORT(qvm);
    REPORT(itc);
    REPORT(btt);
    REPORT(tsq);
    REPORT(indiv_ll);
  }
  return dll;
}

