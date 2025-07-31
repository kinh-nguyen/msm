#include <TMB.hpp>
#include "ktools.hpp"

#define N_PAR 6
#define N_Q 6
#define eps CppAD::numeric_limits<double>::epsilon()

using Eigen::seqN;
using Eigen::Triplet;

template <class Type>
struct makeQ {
  Eigen::SparseMatrix<Type> operator()(const vector<Type> &q_rs, int size) {
    Eigen::SparseMatrix<Type> Q(size, size);
    vector<Type> row_sums = vector<Type>::Zero(size);
    auto add = [&](int r, int c, Type val) {
      if (r < size && c < size) {
        Q.coeffRef(r, c) += val;
        row_sums(r) += val;
      }
    };
    add(0, 1, q_rs(0));
    add(0, 2, q_rs(1));
    add(1, 2, q_rs(2));
    add(2, 3, q_rs(3));
    add(2, 4, q_rs(4));
    add(3, 5, q_rs(5));
    add(4, 5, q_rs(5));
    add(5, 3, q_rs(3));
    add(5, 4, q_rs(4));
    for (int i = 0; i < size; ++i) Q.coeffRef(i, i) = -row_sums(i);
    Q.makeCompressed();
    return Q;
  }
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
  DATA_IVECTOR(dur);
  DATA_VECTOR(n);
  DATA_IVECTOR(fit);
  DATA_VECTOR(afs);
  DATA_VECTOR(afm);
  DATA_VECTOR(aai); // centered

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

  PARAMETER_VECTOR(gp_b);
  dll -= dnorm(gp_b, Type(0), Type(1), true).sum();

  vector<Type> qrs(N_PAR), eta(N_PAR);
  
  auto fill_qrs = [&](Type age, int i) 
  {
    for (int z = 0; z < 3; z++) {
      Type tau_tmp = exp(log_tau[z] + tau_aai[z] * aai[i]);
      eta[z] = itc[z] + logSHASHz(age, mu[z], sigma[z], nu[z], tau_tmp);
    }
    for (int z = 3; z < N_PAR; z++) {
      eta[z] = itc[z] + gp_b[z - 3] * age;
    }
    Type tx = log(1 + exp(age - afs[i]));
    Type tm = log(1 + exp(age - afm[i]));
    eta[2] += b_tx * tx;
    eta({3, 4, 5}) += b_tm * tm;
    qrs = exp(eta);
  };
  
  makeQ<Type> aQ;
  sparse_matrix_exponential::config<Type> cfg_me;
  cfg_me.trace = false;

  for (int i = 0; i < A.size(); i++)
  {
    int m_size = fit[i];
    vector<Type> vp(m_size); vp.setZero();
    vp(A[i]) = 1.;
    for (int j = 0; j < dur[i]; j++) {
      fill_qrs(start[i] + Type(j), i);
      Eigen::SparseMatrix<Type> Q_j = aQ(qrs, m_size);
      sparse_matrix_exponential::expm_generator<Type> p_gen(Q_j, cfg_me);
      vp = p_gen(vp);
    }
    dll -= n[i] * log(vp(Z[i]) + eps);
  }

  SIMULATE
  {
    int n_sim = 80;
    array<Type> PP(N_Q, N_Q, n_sim);
    PP.setZero();
    for (int i = 0; i < n_sim; i++) {
      for (int z = 0; z < 3; z++)
        eta[z] = itc[z] + logSHASHz(Type(i), mu[z], sigma[z], nu[z], tau[z]);
      for (int z = 3; z < N_PAR; z++)
        eta[z] = itc[z] + gp_b[z - 3] * Type(i);
      eta[2] += b_tx;          // 1 year since x
      eta({3, 4, 5}) += b_tm ; // 1 year since m
      qrs = exp(eta);
      matrix<Type> Qd = aQ(qrs, 6).toDense();
      PP.col(i) = atomic::expm(Qd);
    }
    REPORT(PP);
  }
  return dll;
}

