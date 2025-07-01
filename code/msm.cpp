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
struct Qmatrix {
    matrix<Type> Q = matrix<Type>(N_Q, N_Q);
    Qmatrix () {};
    matrix<Type> operator() (vector<Type> q_rs, int size) {
      Q.setZero();
      Q(0, 1) = q_rs(0); // debut
      Q(0, 2) = q_rs(1); // marriage from virgin
      Q(1, 2) = q_rs(2); // marriage from debut
      Q(2, {3, 4, 5}) = q_rs({3, 4, 5}); 
      Q({3, 4, 5}, 6) = q_rs({6, 6, 6}); 
      Q(6, {3, 4, 5}) = q_rs({3, 4, 5});
      Q.diagonal() = Type(-1) * Q.rowwise().sum();
      matrix<Type> q_sub = Q.block(0, 0, size, size);
      return q_sub;
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

  Qmatrix<Type> Q;
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
    Type a = dstart[i];
    Type b = end[i];
    Type Delta = b - a;
    int M = CppAD::Integer(Delta);
    int m_size = fit[i];
    matrix<Type> Qk(m_size, m_size);
    matrix<Type> Rbar(m_size, m_size);
    Rbar.setZero();
    Type lambda_max = 0;
    for (int m = 0; m < M; m++) {
      Type t_m = dstart[i] + (Delta / Type(M)) * (Type(m) + Type(0.5));
      fill_qrs(t_m, i);
      Qk = Q(qrs, m_size);
      Type qmax = Qk.diagonal().cwiseAbs().maxCoeff();
      lambda_max = 0.5 * (lambda_max + qmax + CppAD::abs(lambda_max - qmax));
      matrix<Type> Rk = Qk / lambda_max;
      for (int e = 0; e < m_size; m++) Rk(e, e) += 1;
      Rbar = Rbar + Rk;
    }
    Rbar /= Type(M);
    // sum series P = sum_{k=0}^Kmax e^{-λΔ}(λΔ)^k/k! * Rbar^k
    int Kmax = M + 1;
    matrix<Type> term(m_size, m_size),
        Psum(m_size, m_size);
    term.setIdentity();
    Psum.setZero();
    Type weight = exp(-lambda_max * Delta);
    for (int k = 0; k <= Kmax; k++) {
      Psum += weight * term;
      weight *= (lambda_max * Delta) / Type(k + 1);
      term = Rbar * term;
    }
    dll -= log(Psum(A[i], Z[i]) + Type(1e-16));
  }

  SIMULATE
  {
    // int n_sim = sim_data.rows();
    // array<Type> PP(N_Q, N_Q, n_sim); PP.setZero();
    // for (int i = 0; i < n_sim; i++) {
    //   for (int z = 0; z < N_PAR; z++)
    //     eta[z] = itc[z] + logSHASHz(sim_data(i, 0), mu[z], sigma[z], nu[z], tau[z]);
    //   eta[2] += b_tx * sim_data(i, 1);
    //   eta({3, 4, 5, 6}) += b_tm * sim_data(i, 2);
    //   qrs = exp(eta);
    //   PP.col(i) = KM(qrs, 7);
    // }
    // REPORT(PP);
  }
  return dll;
}

