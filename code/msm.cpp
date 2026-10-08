#include <TMB.hpp>

// Eight-state model (see DECISIONS.md)
//   0 V   never had sex, never married
//   1 X   had sex, never married
//   2 M   in first union
//   3 D1  first union ended by divorce/separation, not remarried
//   4 W1  first union ended by widowhood, not remarried
//   5 R   in a union of order >= 2
//   6 D2  a union of order >= 2 ended by divorce/separation
//   7 W2  a union of order >= 2 ended by widowhood
// Six distinct rates:
//   q0 V->X debut, q1 V->M marriage at debut, q2 X->M marriage after debut,
//   q3 divorce, q4 widowhood, q5 remarriage.
// For the HIV model: U = M + R, D = D1 + D2, W = W1 + W2 (exact under the ties).
#define N_PAR 6
#define N_Q 8
#define eps CppAD::numeric_limits<double>::epsilon()

using Eigen::seqN;
using Eigen::Triplet;

template <class Type>
struct makeQ {
  // size truncates the chain to states 0..size-1; with the state order above
  // every observed cell is reached only through lower-indexed states, so
  // size 3 (pre-marriage), 5 (end in M, D1, W1) and 8 (end in R, D2, W2)
  // give the exact cell probabilities.
  Eigen::SparseMatrix<Type> operator()(const vector<Type> &q, int size) {
    Eigen::SparseMatrix<Type> Q(size, size);
    vector<Type> row_sums = vector<Type>::Zero(size);
    auto add = [&](int r, int c, Type val) {
      if (r < size && c < size) {
        Q.coeffRef(r, c) += val;
        row_sums(r) += val;
      }
    };
    add(0, 1, q(0)); // V  -> X   debut
    add(0, 2, q(1)); // V  -> M   marriage at debut
    add(1, 2, q(2)); // X  -> M   marriage after debut
    add(2, 3, q(3)); // M  -> D1  divorce, first union
    add(2, 4, q(4)); // M  -> W1  widowhood, first union
    add(3, 5, q(5)); // D1 -> R   remarriage
    add(4, 5, q(5)); // W1 -> R   tie: = D1 -> R
    add(5, 6, q(3)); // R  -> D2  tie: = M -> D1
    add(5, 7, q(4)); // R  -> W2  tie: = M -> W1
    add(6, 5, q(5)); // D2 -> R   tie: 3rd+ unions unobserved
    add(7, 5, q(5)); // W2 -> R   tie
    for (int i = 0; i < size; ++i) Q.coeffRef(i, i) = -row_sums(i);
    Q.makeCompressed();
    return Q;
  }
};

template <class Type>
Type logSHASHz(Type t, Type mu, Type sigma, Type nu, Type tau)
{
  // log hazard of a sinh-arcsinh distribution on log age
  Type x = log(t),
       z = (x - mu) / (sigma * tau),
       tau_asinh_nu = tau * log(z + sqrt(z * z + 1)) - nu,
       c = cosh(tau_asinh_nu),
       r = sinh(tau_asinh_nu);
  Type logres = -log(sigma) - 0.5 * log(2 * M_PI) - 0.5 * log(1 + (z * z)) + log(c) - 0.5 * (r * r) - x;
  Type logp = log(1.0 - pnorm(r) + eps);
  return logres - logp;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  parallel_accumulator<Type> dll(this);

  // data, one row per aggregated episode
  DATA_IVECTOR(A);      // start state
  DATA_IVECTOR(Z);      // end state (observed cell)
  DATA_VECTOR(start);   // age at start of episode
  DATA_IVECTOR(dur);    // years in episode
  DATA_VECTOR(n);       // Kish weight divided by the design effect
  DATA_IVECTOR(fit);    // chain size for this episode: 3, 5 or 8
  DATA_VECTOR(afs);     // age at first sex (100 if none)
  DATA_VECTOR(afm);     // age at first marriage (100 if none)
  DATA_VECTOR(coh);     // birth cohort, centred, in decades

  DATA_SCALAR(age_c);   // centring of age in the post-marriage rates
  DATA_SCALAR(tm_c);    // centring of duration since first marriage
  DATA_SCALAR(itc_mu);  // prior mean of the post-marriage intercepts
  DATA_SCALAR(itc_sd);  // prior sd of the post-marriage intercepts
  DATA_SCALAR(afs_ref); // reference person for SIMULATE
  DATA_SCALAR(afm_ref);
  DATA_SCALAR(coh_ref);

  // pre-marriage rates: SHASH hazard in age
  PARAMETER_VECTOR(mu);        // location (log age)
  PARAMETER_VECTOR(log_sigma); // spread
  PARAMETER_VECTOR(log_nu);    // skewness
  PARAMETER_VECTOR(log_tau);   // tail weight
  PARAMETER_VECTOR(b_coh);     // linear cohort shift of mu, per decade
  PARAMETER(b_tx);             // X -> M: time since debut
  vector<Type> sigma = exp(log_sigma), nu = exp(log_nu), tau = exp(log_tau);

  dll -= dnorm(mu, Type(3), Type(1), true).sum(); // log(20)
  dll -= dnorm(nu, Type(0), Type(1), true).sum() + log_nu.sum();
  dll -= dnorm(sigma, Type(0), Type(1), true).sum() + log_sigma.sum();
  dll -= dnorm(tau, Type(0), Type(1), true).sum() + log_tau.sum();
  dll -= dnorm(b_coh, Type(0), Type(1), true).sum();
  dll -= dnorm(b_tx, Type(0), Type(1), true);

  // intercepts; the first three are mapped to 0 (the SHASH carries the level)
  PARAMETER_VECTOR(itc);
  for (int z = 3; z < N_PAR; z++) dll -= dnorm(itc[z], itc_mu, itc_sd, true);

  // post-marriage rates: log-linear in age, duration since first marriage, cohort
  PARAMETER_VECTOR(gp_b);       // age slope, per year
  PARAMETER_VECTOR(b_tm);       // duration since first marriage, per year
  PARAMETER_VECTOR(b_coh_post); // cohort, per decade (mapped off in single-survey countries)
  dll -= dnorm(gp_b, Type(0), Type(1), true).sum();
  dll -= dnorm(b_tm, Type(0), Type(1), true).sum();
  dll -= dnorm(b_coh_post, Type(0), Type(1), true).sum();

  vector<Type> qrs(N_PAR), eta(N_PAR);

  auto fill_qrs = [&](Type age, Type afs_i, Type afm_i, Type coh_i) {
    Type tm = log(1 + exp(age - afm_i)); // softplus: duration since first marriage
    Type tx = log(1 + exp(age - afs_i)); // softplus: duration since debut
    for (int z = 0; z < 3; z++)
      eta[z] = itc[z] + logSHASHz(age, mu[z] + b_coh[z] * coh_i, sigma[z], nu[z], tau[z]);
    eta[2] += b_tx * tx;
    for (int z = 3; z < N_PAR; z++)
      eta[z] = itc[z] + gp_b[z - 3] * (age - age_c) + b_tm[z - 3] * (tm - tm_c)
             + b_coh_post[z - 3] * coh_i;
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
      fill_qrs(start[i] + Type(j), afs[i], afm[i], coh[i]);
      Eigen::SparseMatrix<Type> Q_j = aQ(qrs, m_size);
      sparse_matrix_exponential::expm_generator<Type> p_gen(Q_j, cfg_me);
      vp = p_gen(vp);
    }
    dll -= n[i] * log(vp(Z[i]) + eps);
  }

  SIMULATE
  {
    // annual transition matrices and rates by age for one reference person
    int n_sim = 80;
    array<Type> PP(N_Q, N_Q, n_sim);
    matrix<Type> QQ(n_sim, N_PAR);
    PP.setZero(); QQ.setZero();
    for (int a = 1; a < n_sim; a++) {
      fill_qrs(Type(a), afs_ref, afm_ref, coh_ref);
      for (int z = 0; z < N_PAR; z++) QQ(a, z) = qrs[z];
      matrix<Type> Qd = aQ(qrs, N_Q).toDense();
      PP.col(a) = atomic::expm(Qd);
    }
    REPORT(PP);
    REPORT(QQ);
  }
  return dll;
}
