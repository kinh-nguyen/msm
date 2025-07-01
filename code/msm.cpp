#include <TMB.hpp>
#include "ktools.hpp"

#define N_PAR 7
#define N_Q 7
#define eps CppAD::numeric_limits<double>::epsilon()

using Eigen::seqN;

template <class T> 
struct PQ {
    matrix<T> Q = matrix<T>(N_Q, N_Q);
    PQ () {};
    matrix<T> operator() (vector<T> q_rs, int size) {
      Q.setZero();
      Q(0, 1) = q_rs(0); // debut
      Q(0, 2) = q_rs(1); // marriage from virgin
      Q(1, 2) = q_rs(2); // marriage from debut
      Q(2, {3, 4, 5}) = q_rs({3, 4, 5}); 
      Q({3, 4, 5}, 6) = q_rs({6, 6, 6}); 
      Q(6, {3, 4, 5}) = q_rs({3, 4, 5});
      Q.diagonal() = T(-1) * Q.rowwise().sum();
      matrix<T> q_sub = Q.block(0, 0, size, size);
      return atomic::expm(q_sub);
    };
};

template <class Type>
Type logSHASHz(Type t, Type mu, Type sigma, Type nu, Type tau)
{
  Type x = log(t),
       z = (x - mu) / (sigma * tau),
       tau_asinh_nu = tau * log(z + sqrt(z * z + 1)) - nu,
       c = cosh(tau_asinh_nu),
       r = sinh(tau_asinh_nu);
  Type logres = -log(sigma) - 0.5 * log(2 * M_PI) - 0.5 * log(1 + (z * z)) + log(c) - 0.5 * (r * r) - x;
  Type logp = log(1.0 - pnorm(r) + eps);
  Type log_hz = logres - logp;
  return log_hz;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
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
  DATA_VECTOR(aai); // centered
  DATA_VECTOR(tm);
  DATA_VECTOR(tx);

  DATA_MATRIX(sim_data);

  PARAMETER_VECTOR(mu);
  PARAMETER_VECTOR(log_sigma);
  vector<Type> sigma = exp(log_sigma);
  PARAMETER_VECTOR(log_nu); // no need to be +, but an informative prior to have nu > 0
  vector<Type> nu = exp(log_nu);
  PARAMETER_VECTOR(log_tau); // now the intercept
  vector<Type> tau = exp(log_tau);
  PARAMETER_VECTOR(tau_aai); // one for each transition, small deviation
  dll -= dnorm(tau_aai, Type(0), Type(1), true).sum();
  
  dll -= dnorm(mu, Type(3), Type(1), true).sum(); // log(20)
  dll -= dnorm(nu, Type(0), Type(1), true).sum() + log_nu.sum();
  dll -= dnorm(sigma, Type(0), Type(1), true).sum() + log_sigma.sum();
  dll -= dnorm(tau, Type(0), Type(1), true).sum() + log_tau.sum();

  PARAMETER(b_tx);
  dll -= dnorm(b_tx, Type(0), Type(1), true);

  PARAMETER(b_tm);
  dll -= dnorm(b_tm, Type(0), Type(1), true);

  PARAMETER_VECTOR(itc);
  dll -= dnorm(itc, Type(0), Type(1), true).sum();

  PQ<Type> KM;
  vector<Type> qrs(N_PAR), eta(N_PAR);

  auto fill_qrs = [&](Type age, int i) 
  {
    for (int z = 0; z < N_PAR; z++) {
      Type tau_tmp = exp(log_tau[z] + tau_aai[z] * aai[i]);
      eta[z] = itc[z] + logSHASHz(age, mu[z], sigma[z], nu[z], tau_tmp);
    }
    eta[2] += b_tx * tx[i];
    eta({3, 4, 5, 6}) += b_tm * tm[i];
    qrs = exp(eta);
  };

  for (int i = 0; i < A.size(); i++) 
  {
    int m_size = fit[i];
    if (m_size == 3) {
      fill_qrs(dstart[i], i);
      dll -= n[i] * log(KM(qrs, m_size)(A[i], Z[i]) + eps);
    } else {
      vector<Type> vp(m_size);
      vp.setZero();
      vp(A[i]) = 1.;
      for (int j = istart[i]; j < end[i]; j++)
      {
        fill_qrs(Type(j), i);
        matrix<Type> km = KM(qrs, m_size);
        vp = vp.matrix() * km;
      }
      dll -= n[i] * log(vp(Z[i]) + eps);
    }
  }

  SIMULATE
  {
    int n_sim = sim_data.rows();
    array<Type> PP(N_Q, N_Q, n_sim); PP.setZero();
    for (int i = 0; i < n_sim; i++) {
      for (int z = 0; z < N_PAR; z++)
        eta[z] = itc[z] + logSHASHz(sim_data(i, 0), mu[z], sigma[z], nu[z], tau[z]);
      eta[2] += b_tx * sim_data(i, 1);
      eta({3, 4, 5, 6}) += b_tm * sim_data(i, 2);
      qrs = exp(eta);
      PP.col(i) = KM(qrs, 7);
    }
    REPORT(PP);
  }
  return dll;
}

