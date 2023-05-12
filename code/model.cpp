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


  // base rate | intercept
  PARAMETER_VECTOR(intercepts);
  prior -= dnorm(intercepts, prior_base(0), prior_base(1), true).sum();

  // age's coeff | time in the hazard
  PARAMETER_VECTOR(beta_t);
  prior -= dnorm(beta_t, prior_t(0), prior_t(1), true).sum();

  // cc's coeff | random intercept
  PARAMETER_VECTOR(cc0);
  PARAMETER_VECTOR(cc1);
  PARAMETER_VECTOR(cc2);
  PARAMETER_VECTOR(cc3);
  PARAMETER_VECTOR(cc4);
  PARAMETER_VECTOR(cc5);
  PARAMETER_VECTOR(cc6);
  prior -= dnorm(cc0, prior_cc(0), prior_cc(1), true).sum();
  prior -= dnorm(cc1, prior_cc(0), prior_cc(1), true).sum();
  prior -= dnorm(cc2, prior_cc(0), prior_cc(1), true).sum();
  prior -= dnorm(cc3, prior_cc(0), prior_cc(1), true).sum();
  prior -= dnorm(cc4, prior_cc(0), prior_cc(1), true).sum();
  prior -= dnorm(cc5, prior_cc(0), prior_cc(1), true).sum();
  prior -= dnorm(cc6, prior_cc(0), prior_cc(1), true).sum();
  matrix<Type> ccmat(N_CC, N_PAR);
  ccmat.col(0) = cc0;
  ccmat.col(1) = cc1;
  ccmat.col(2) = cc2;
  ccmat.col(3) = cc3;
  ccmat.col(4) = cc4;
  ccmat.col(5) = cc5;
  ccmat.col(6) = cc6;

  PQ<Type> KM;

  vector<Type> ll(A.size());
  vector<Type> eta(N_PAR);
  for (int i = 0; i < A.size(); i++) {
    for (int j = 0; j < N_PAR; j++)
      eta[j] = exp(intercepts[j] + beta_t[j]*start[i] + ccmat(cid[i], j));
    ll[i] = n[i] * KM(eta, start[i], end[i])(A[i], Z[i]);
  }
  REPORT(prior);
  REPORT(ll);
  dll += prior - ll.sum();
  REPORT(intercepts);
  REPORT(beta_t);
  REPORT(cc0);
  REPORT(cc1);
  REPORT(cc2);
  REPORT(cc3);
  REPORT(cc4);
  REPORT(cc5);
  REPORT(cc6);
  return dll;
}
