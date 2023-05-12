#include <TMB.hpp>

#define N_PAR 7
#define N_Q 7
#define N_CC 37

using Eigen::seqN;

template <class T> 
struct PQ {
    matrix<T> 
        qM = matrix<T>(N_Q, N_Q),
        pM = matrix<T>(N_Q, N_Q);
    PQ () {};
    matrix<T> operator() (vector<T> eta, T start, T end, bool isLog = true) {
      qM.setZero(); pM.setZero();
      qM(0, 1) = eta(0); // debut
      qM(0, 2) = eta(1); // marriage from virgin
      qM(1, 2) = eta(2); // marriage from debut
      qM(2, {3,4,5}) = eta({3,4,5}); // marriage dissolution
      qM({3,4,5}, 6) = eta({6,6,6}); // disso > remarried
      qM(6, {3,4,5}) = eta({3,4,5}); // remarried > disso = married > disso, we could add a(three) scaling parameter as well?
      qM.diagonal() = T(-1) * qM.rowwise().sum();
      qM *= (end - start);
      pM = expm(qM);
      if (isLog) pM = log(pM.array());
      return pM;
    };
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);
  Type prior = 0.0;

  // data
  DATA_IVECTOR(A);
  DATA_IVECTOR(Z);
  DATA_IVECTOR(start);
  DATA_IVECTOR(end);
  DATA_VECTOR(n);
  DATA_IVECTOR(cid);
  
  // priors
  DATA_VECTOR(prior_base);
  DATA_VECTOR(prior_t);
  DATA_VECTOR(prior_cc);

  // Coefs
  PARAMETER_VECTOR(betas);

  // base rate | intercept
  vector<Type> intercepts = betas(seqN(0, N_PAR));
  prior -= dnorm(intercepts, prior_base(0), prior_base(1), true).sum();

  // age's coeff | time in the hazard
  vector<Type> beta_t = betas(seqN(N_PAR, N_PAR));
  prior -= dnorm(beta_t, prior_t(0), prior_t(1), true).sum();

  // cc's coeff | random intercept
  matrix<Type> ccmat = betas(seqN(N_PAR*2, N_CC*N_PAR)).reshaped(N_CC, N_PAR);
  for (int i = 0; i < N_PAR; i++) {
    vector<Type> cci = ccmat.col(i);
    prior -= dnorm(cci, prior_cc(0), prior_cc(1), true).sum();
  }

  PQ<Type> KM;

  vector<Type> eta(N_PAR);
  for (int i = 0; i < A.size(); i++) {
    for (int j = 0; j < N_PAR; j++)
      eta[j] = exp(intercepts[j] + beta_t[j]*start[i] + ccmat(cid[i], j));
    dll -= n[i] * KM(eta, start[i], end[i])(A[i], Z[i]);
  }
  dll += prior;

  REPORT(betas);
  return dll;
}
