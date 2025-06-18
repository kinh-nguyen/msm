#include <TMB.hpp>
#include "ktools.hpp"

#define N_PAR 3
#define N_Q 3
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
  
  DATA_VECTOR(tx);

  DATA_MATRIX(sim_data);

  DATA_MATRIX(Bspline);

  PARAMETER_VECTOR(VX);                      
  PARAMETER_VECTOR(VM);
  PARAMETER_VECTOR(XM);

  vector<Type>
      VXv = Bspline * VX,
      VMv = Bspline * VM,
      XMv = Bspline * XM
      ;

  PARAMETER_VECTOR(log_k); // smooth penalty
  vector<Type> k(exp(log_k));
  dll -= dnorm(k, Type(0), Type(1), true).sum() + log_k.sum();

  DATA_SPARSE_MATRIX(penalty);
  SparseMatrix<Type>
      Q0 = k[0] * penalty,
      Q1 = k[1] * penalty,
      Q2 = k[2] * penalty
      ;

  dll += GMRF(Q0)(VX) +
         GMRF(Q1)(VM) +
         GMRF(Q2)(XM)
         ;

  PARAMETER(b_tx);
  dll -= dnorm(b_tx, Type(0), Type(1), true);

  PARAMETER_VECTOR(itc);
  dll -= dnorm(itc, Type(0), Type(1), true).sum();

  PQ<Type> KM;
  vector<Type> qrs(N_PAR);

  vector<Type> indiv_ll(A.size());
  indiv_ll.setZero();

  auto fill_qrs = [&](int age, int i)
  {
    qrs[0] = exp(itc[0] + VXv[age]); 
    qrs[1] = exp(itc[1] + VMv[age]); 
    qrs[2] = exp(itc[2] + XMv[age] + b_tx * tx[i]);
  };

  matrix<Type>
      Pm(N_Q, N_Q), Qm(N_Q, N_Q);
  Type ll_val = 0;

  for (int i = 0; i < A.size(); i++)
  {
    if (fit[i] == 0)
    {
      fill_qrs(istart[i], i);
      Qm = KM(qrs, false, false);
      ll_val = log(Qm(A[i], Z[i]) + eps);
    }
    else if (fit[i] == 1)
    {
      fill_qrs(istart[i], i);
      Pm = KM(qrs, false, false);
      ll_val = log(Pm(A[i], Z[i]) + eps);
    }
    dll -= n[i] * ll_val;
    indiv_ll[i] = -n[i] * ll_val;
  }

  SIMULATE
  {
    int n_sim = sim_data.rows();
    array<Type> PP(N_Q, N_Q, n_sim); PP.setZero();
    for (int i = 0; i < n_sim; i++) {
      int age = CppAD::Integer(sim_data(i, 0));
      qrs[0] = exp(itc[0] + VXv[age]);
      qrs[1] = exp(itc[1] + VMv[age]);
      qrs[2] = exp(itc[2] + XMv[age] + b_tx * tx[i]);
      Pm = KM(qrs, false, false);
      PP.col(i) = Pm;
    }
    REPORT(PP);
  }
  return dll;
}

