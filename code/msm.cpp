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
  DATA_VECTOR(tm);

  DATA_MATRIX(sim_data);

  DATA_MATRIX(Bspline);

  PARAMETER_VECTOR(VX);                      
  PARAMETER_VECTOR(VM);
  PARAMETER_VECTOR(XM);
  PARAMETER_VECTOR(MS);
  PARAMETER_VECTOR(MD);
  PARAMETER_VECTOR(MW);
  PARAMETER_VECTOR(UR);

  vector<Type>
      VXv = Bspline * VX,
      VMv = Bspline * VM,
      XMv = Bspline * XM,
      MSv = Bspline * MS,
      MDv = Bspline * MD,
      MWv = Bspline * MW, 
      URv = Bspline * UR;

  PARAMETER_VECTOR(log_k); // smooth penalty
  vector<Type> k(exp(log_k));
  dll -= dnorm(k, Type(0), Type(1), true).sum() + log_k.sum();

  DATA_SPARSE_MATRIX(penalty);
  SparseMatrix<Type>
      Q0 = k[0] * penalty,
      Q1 = k[1] * penalty,
      Q2 = k[2] * penalty,
      Q3 = k[3] * penalty,
      Q4 = k[4] * penalty,
      Q5 = k[5] * penalty,
      Q6 = k[6] * penalty;

  dll += GMRF(Q0)(VX) +
         GMRF(Q1)(VM) +
         GMRF(Q2)(XM) +
         GMRF(Q3)(MS) +
         GMRF(Q4)(MD) +
         GMRF(Q5)(MW) +
         GMRF(Q6)(UR);

  PARAMETER(b_tx);
  PARAMETER(b_tm);
  dll -= dnorm(b_tx, Type(0), Type(1), true);
  dll -= dnorm(b_tm, Type(0), Type(1), true);

  PARAMETER_VECTOR(itc);
  dll -= dnorm(itc, Type(0), Type(1), true).sum();

  PARAMETER(logit_pi);
  Type pi = invlogit(logit_pi);
  dll -= dbeta(pi, Type(5), Type(2), true);
  dll -= log(pi) + log(1 - pi);

  PQ<Type> KM;
  vector<Type> qrs(N_PAR), srs(N_PAR);

  vector<Type> indiv_ll(A.size());
  indiv_ll.setZero();

  auto fill_qrs = [&](int age, int i)
  {
    qrs[0] = exp(itc[0] + VXv[age]); // pop 1, force pop 2 does not have VX
    qrs[1] = Type(0); // pop 2, force pop 1 does not have VM
    qrs[2] = exp(itc[2] + XMv[age] + b_tx * tx[i]);
    qrs[3] = exp(itc[3] + MSv[age] + b_tm * tm[i]);
    qrs[4] = exp(itc[4] + MDv[age] + b_tm * tm[i]);
    qrs[5] = exp(itc[5] + MWv[age] + b_tm * tm[i]);
    qrs[6] = exp(itc[6] + URv[age]);
    srs[0] = Type(0); // force pop 2 does not have VX
    srs[1] = exp(itc[1] + VMv[age]); // pop 2, force pop 1 does not have VM
    srs[2] = Type(0);
    srs[3] = exp(itc[3] + MSv[age] + b_tm * tm[i]);
    srs[4] = exp(itc[4] + MDv[age] + b_tm * tm[i]);
    srs[5] = exp(itc[5] + MWv[age] + b_tm * tm[i]);
    srs[6] = exp(itc[6] + URv[age]);
  };

  matrix<Type> 
    Pm(N_Q, N_Q), Qm(N_Q, N_Q), cumPm(N_Q, N_Q),
    Pm1(N_Q, N_Q), Qm1(N_Q, N_Q), cumPm1(N_Q, N_Q);
  Type ll_val = 0;

  for (int i = 0; i < A.size(); i++)
  {
    if (fit[i] == 0)
    {
      fill_qrs(istart[i], i);
      Qm = KM(qrs, false, false);
      Qm1 = KM(srs, false, false);
      ll_val = log(pi * Qm(A[i], Z[i]) + (1 - pi) * Qm1(A[i], Z[i]) + eps);
    }
    else if (fit[i] == 1)
    {
      fill_qrs(istart[i], i);
      Pm = KM(qrs, false, false);
      Pm1 = KM(srs, false, false);
      ll_val = log(pi * Pm(A[i], Z[i]) + (1 - pi) * Pm1(A[i], Z[i]) + eps);
    }
    else if (fit[i] == 2)
    { // no need two pops here at there are no VM/VX
      cumPm = matrix<Type>::Identity(N_Q, N_Q);
      cumPm1 = matrix<Type>::Identity(N_Q, N_Q);
      for (int t = istart[i]; t < end[i]; t++)
      {
        fill_qrs(t, i);
        Pm = KM(qrs, false, false);
        Pm1 = KM(srs, false, false);
        cumPm = cumPm * Pm;
        cumPm1 = cumPm1 * Pm1;
      }
      ll_val = log(pi * cumPm(A[i], Z[i]) + (1 - pi) * cumPm1(A[i], Z[i]) + eps);
    }
    dll -= n[i] * ll_val;
    indiv_ll[i] = -n[i] * ll_val;
  }

  SIMULATE
  {
    int n_sim = sim_data.rows();
    array<Type> PP(N_Q, N_Q, n_sim); PP.setZero();
    array<Type> PP1(N_Q, N_Q), PP2(N_Q, N_Q); 
    PP1.setZero();
    PP2.setZero();
    for (int i = 0; i < n_sim; i++) {
      int age = CppAD::Integer(sim_data(i, 0));
      qrs[0] = exp(itc[0] + VXv[age]);
      qrs[1] = Type(0);
      qrs[2] = exp(itc[2] + XMv[age] + b_tx * sim_data(i, 1));
      qrs[3] = exp(itc[3] + MSv[age] + b_tm * sim_data(i, 2));
      qrs[4] = exp(itc[4] + MDv[age] + b_tm * sim_data(i, 2));
      qrs[5] = exp(itc[5] + MWv[age] + b_tm * sim_data(i, 2));
      qrs[6] = exp(itc[6] + URv[age]);
      srs[0] = Type(0);
      srs[1] = exp(itc[1] + VMv[age]);
      srs[2] = Type(0);
      srs[3] = exp(itc[3] + MSv[age] + b_tm * sim_data(i, 2));
      srs[4] = exp(itc[4] + MDv[age] + b_tm * sim_data(i, 2));
      srs[5] = exp(itc[5] + MWv[age] + b_tm * sim_data(i, 2));
      srs[6] = exp(itc[6] + URv[age]);
      PP1 = KM(qrs, false, false);
      PP2 = KM(srs, false, false);
      PP.col(i) = pi * PP1 + (1 - pi) * PP2;
    }
    REPORT(PP);
  }
  return dll;
}

