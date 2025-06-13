#include <TMB.hpp>
#include "ktools.hpp"

#define N_PAR 7
#define N_Q 7
#define eps CppAD::numeric_limits<double>::epsilon()

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

template <class Type>
Type logSHASHz(Type t, Type mu, Type sigma, Type nu, Type tau)
{
  Type x = log(t),
       z = (x - mu) / sigma,
       tau_asinh_nu = tau * log(z + sqrt(z * z + 1)) - nu,
       c = cosh(tau_asinh_nu),
       r = sinh(tau_asinh_nu);
  // SHASHo2 with sigma' = sigma . tau
  Type logres = -log(sigma) - 0.5 * log(2 * M_PI) - 0.5 * log(1 + (z * z)) + log(c) - 0.5 * (r * r) - x;
  Type logp = log(1.0 - pnorm(r) + eps);
  Type log_hz = logres - logp;
  return log_hz;
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
  DATA_IVECTOR(end);
  DATA_VECTOR(n);
  DATA_IVECTOR(fit);
  DATA_INTEGER(n_age);
  DATA_VECTOR(tx);
  DATA_VECTOR(tm);

  DATA_MATRIX(sim_data);

  PARAMETER_VECTOR(mu);
  DATA_VECTOR(prior_mu);
  vector<Type> exp_mu = exp(mu);
  dll -= dnorm(mu, prior_mu(0), prior_mu(1), true).sum();

  PARAMETER_VECTOR(log_sigma);
  DATA_VECTOR(prior_sigma);
  vector<Type> sigma = exp(log_sigma);
  dll -= dnorm(sigma, prior_sigma(0), prior_sigma(1), true).sum() + log_sigma.sum();

  PARAMETER_VECTOR(nu);
  DATA_VECTOR(prior_nu);
  dll -= dnorm(nu, prior_nu(0), prior_nu(1), true).sum();
  
  PARAMETER_VECTOR(log_tau);
  DATA_VECTOR(prior_tau);
  vector<Type> tau = exp(log_tau);
  dll -= dnorm(tau, prior_tau(0), prior_tau(1), true).sum() + log_tau.sum();

  DATA_VECTOR(prior_coef);
  
  PARAMETER_VECTOR(itc);
  dll -= dnorm(itc, prior_coef(0), prior_coef(1), true).sum();

  PARAMETER(b_tx);
  PARAMETER(b_tm);
  dll -= dnorm(b_tx, Type(0), Type(1), true);
  dll -= dnorm(b_tm, Type(0), Type(1), true);

  PQ<Type> KM;
  vector<Type> qrs(N_PAR);

  vector<Type> indiv_ll(A.size());
  indiv_ll.setZero();

  auto fill_qrs = [&](Type age, int i)
  {
    for (int p = 0; p < N_PAR; p++)
    {
      Type lhz = logSHASHz(age, mu[p], sigma[p], nu[p], tau[p]);
      qrs[p] = exp(itc[p] + lhz);
      if (p == 2) 
        qrs[p] *= exp(b_tx * tx[i]); // time since debuted
      if (p > 2)
        qrs[p] *= exp(b_tm * tm[i]); // time since married
    }
  };

  matrix<Type> Pm(N_Q, N_Q), Qm(N_Q, N_Q), cumPm(N_Q, N_Q);
  Type ll_val = 0;

  for (int i = 0; i < A.size(); i++)
  {
    if (fit[i] == 0)
    {
      fill_qrs(dstart[i], i);
      Qm = KM(qrs, false, true);
      ll_val = log(Qm(A[i], Z[i]) + eps);
    }
    else if (fit[i] == 1)
    {
      fill_qrs(dstart[i], i);
      Pm = KM(qrs, false, false);
      ll_val = log(Pm(A[i], Z[i]) + eps);
    }
    else if (fit[i] == 2)
    {
      cumPm = matrix<Type>::Identity(N_Q, N_Q);
      for (int t = istart[i]; t < end[i]; t++)
      {
        fill_qrs(Type(t), i);
        Pm = KM(qrs, false, false);
        cumPm = cumPm * Pm;
      }
      ll_val = log(cumPm(A[i], Z[i]) + eps);
    }
    dll -= n[i] * ll_val;
    indiv_ll[i] = -n[i] * ll_val;
  }

  SIMULATE
  {
    int n_sim = sim_data.rows();
    array<Type> PP(N_Q, N_Q, n_sim); PP.setZero();
    for (int i = 0; i < n_sim; i++) {
        for (int p = 0; p < N_PAR; p++) {
          Type lhz = logSHASHz(sim_data(i, 0), mu[p], sigma[p], nu[p], tau[p]);
          qrs[p] = exp(itc[p] + lhz);
          if (p == 2)
            qrs[p] *= exp(b_tx * sim_data(i, 1)); // time since debuted
          if (p > 2)
            qrs[p] *= exp(b_tm * sim_data(i, 2)); // time since married
        }
        PP.col(i) = KM(qrs, false, false);
    }
   
    REPORT(PP);
    REPORT(exp_mu);
    REPORT(sigma);
    REPORT(tau);
    REPORT(nu);
    REPORT(itc);
    REPORT(indiv_ll);
  }
  return dll;
}

