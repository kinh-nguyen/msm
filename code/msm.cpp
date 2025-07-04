#include <TMB.hpp>
#include "ktools.hpp"

#define N_PAR 7
#define N_Q 7
#define eps CppAD::numeric_limits<double>::epsilon()

using Eigen::seqN;

template <class T> 
struct Qs {
    matrix<T> Q = matrix<T>(N_Q, N_Q);
    Qs () {};
    matrix<T> operator() (vector<T> q_rs) {
      Q.setZero();
      Q(0, 1) = q_rs(0); // debut
      Q(0, 2) = q_rs(1); // marriage from virgin
      Q(1, 2) = q_rs(2); // marriage from debut
      Q(2, {3, 4, 5}) = q_rs({3, 4, 5}); 
      Q({3, 4, 5}, 6) = q_rs({6, 6, 6}); 
      Q(6, {3, 4, 5}) = q_rs({3, 4, 5});
      Q.diagonal() = T(-1) * Q.rowwise().sum();
      return Q;
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
  DATA_VECTOR(start);
  DATA_IVECTOR(end);
  DATA_VECTOR(n);
  DATA_IVECTOR(fit);
  DATA_VECTOR(afs);
  DATA_VECTOR(afm);
  DATA_VECTOR(aai); // centered

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

  Qs<Type> Q;
  vector<Type> qrs(N_PAR), eta(N_PAR);

  auto fill_qrs = [&](Type age, int i) 
  {
    for (int z = 0; z < N_PAR; z++) {
      Type tau_tmp = exp(log_tau[z] + tau_aai[z] * aai[i]);
      eta[z] = itc[z] + logSHASHz(age, mu[z], sigma[z], nu[z], tau_tmp);
    }
    Type tx = log(1 + exp(age - afs[i]));
    Type tm = log(1 + exp(age - afm[i]));
    eta[2] += b_tx * tx;
    eta({3, 4, 5, 6}) += b_tm * tm;
    qrs = exp(eta);
  };

  for (int i = 0; i < A.size(); i++) 
  {
    Type lli = 0;
    if (A[i] == 0) // fit = 3
    {
      fill_qrs(start[i], i);
      Type q0 = qrs[0] + qrs[1];
      lli = -q0;
      if (Z[i] == 0)
        for (int j = 0; j < end[i]; j++)
        {
          fill_qrs(start[i] + Type(j), i);
          lli += -(qrs[0] + qrs[1]);
        }
      if (Z[i] == 1) {
        if (qrs[1] != q0) lli = log(qrs[0]/(qrs[2] - q0) * (exp(-q0) - exp(-qrs[2])));
        if (qrs[1] == q0) lli = log(qrs[0]) - q0;
      }
      if (Z[i] == 2) {
        Type p00 = exp(-q0);
        // if (qrs[1] != q0)
        Type p01 = qrs[0] / (qrs[2] - q0) * (exp(-q0) - exp(-qrs[2]));
        if (qrs[1] == q0)
          p01 = qrs[0] * p00;
        lli = log(1 - p00 - p01);
      }
    }
    if (A[i] == 1) { // fit == 3
      if (Z[i] == 1) 
        for (int j = 0; j < end[i]; j++) {
          fill_qrs(start[i] + Type(j), i);
          lli += -qrs[2];
        }
      else
      {
        fill_qrs(start[i], i);
        lli = log(1 - exp(-qrs[2]));
      }
    }
    if (A[i] == 2 && fit[i] == 6) { // change 6 to 4 in R
      if (Z[i] == 2)
        for (int j = 0; j < end[i]; j++)
        {
          fill_qrs(start[i] + Type(j), i);
          Type q2 = qrs({3, 4, 5}).sum();
          lli += -q2;
        }
      if (Z[i] != 2) {
        Type p22sofar = 1, p2jTotal = 0;
        for (int j = 0; j < end[i]; j++)
        {
          fill_qrs(start[i] + Type(j), i);
          Type q2 = qrs({3, 4, 5}).sum();
          Type p2jt = qrs[Z[i]] / q2 * (1 - exp(-q2));
          p2jTotal += p22sofar * p2jt;
          p22sofar *= exp(-q2);
        }
        lli = log(p2jTotal);
      }
    }
    if (A[i] == 2 && fit[i] == 7) { // change 7 to 5 in R?
      matrix<Type> vp(1, 5);
      vp << 1, 0, 0, 0, 0;
      for (int j = 0; j < end[i]; j++)
      {
        fill_qrs(start[i] + Type(j), i);
        matrix<Type> km = Q(qrs).block(2, 2, 5, 5);
        vp = vp * atomic::expm(km);
      }
      lli = log(vp(Z[i]) + eps);
    }
    dll -= n[i] * lli;
  }

  SIMULATE
  {
    int n_sim = sim_data.rows();
    array<Type> PP(N_Q, N_Q, n_sim); PP.setZero();
    for (int i = 0; i < n_sim; i++) {
      for (int z = 0; z < N_PAR; z++)
        eta[z] = itc[z] + logSHASHz(sim_data(i, 0), mu[z], sigma[z], nu[z], tau[z]);
      eta[2] += b_tx * sim_data(i, 1); // average time since x
      eta({3, 4, 5, 6}) += b_tm * sim_data(i, 2); // averate time since m
      qrs = exp(eta);
      PP.col(i) = atomic::expm(Q(qrs));
    }
    REPORT(PP);
  }
  return dll;
}

