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

  DATA_MATRIX(Mspline);
  PARAMETER_VECTOR(VX);
  PARAMETER_VECTOR(VM);
  PARAMETER_VECTOR(XM);
  PARAMETER_VECTOR(MS);
  PARAMETER_VECTOR(MD);
  PARAMETER_VECTOR(MW);
  PARAMETER_VECTOR(UR);

  DATA_MATRIX(penalty);

  vector<Type>
      VXw = exp(VX) / exp(VX).sum(),
      VMw = exp(VM) / exp(VM).sum(),
      XMw = exp(XM) / exp(XM).sum(),
      MSw = exp(MS) / exp(MS).sum(),
      MDw = exp(MD) / exp(MD).sum(),
      MWw = exp(MW) / exp(MW).sum(),
      URw = exp(UR) / exp(UR).sum(),
      // splines
      VXv = Mspline * VXw,
      VMv = Mspline * VMw,
      XMv = Mspline * XMw,
      MSv = Mspline * MSw,
      MDv = Mspline * MDw,
      MWv = Mspline * MWw,
      URv = Mspline * URw;

  // quad-form
  Type
      VXp = (VX * (penalty * VX)).sum(),
      VMp = (VM * (penalty * VM)).sum(),
      XMp = (XM * (penalty * XM)).sum(),
      MSp = (MS * (penalty * MS)).sum(),
      MDp = (MD * (penalty * MD)).sum(),
      MWp = (MW * (penalty * MW)).sum(),
      URp = (UR * (penalty * UR)).sum();

  PARAMETER(log_k); // smooth penalty
  Type k(exp(log_k));
  dll -= dnorm(k, Type(0), Type(1), true) + log_k;
  dll += 0.5 * k * (VXp + VMp + XMp + MSp + MDp + MWp + URp);

  PARAMETER(b_tx);
  PARAMETER(b_tm);

  dll -= dnorm(b_tx, Type(0), Type(1), true);
  dll -= dnorm(b_tm, Type(0), Type(1), true);

  PQ<Type> KM;
  vector<Type> qrs(N_PAR);

  vector<Type> indiv_ll(A.size());
  indiv_ll.setZero();

  auto fill_qrs = [&](int age, int i)
  {
    qrs[0] = VXv[age];
    qrs[1] = VMv[age];
    qrs[2] = XMv[age] * exp(b_tx * tx[i]);
    qrs[3] = MSv[age] * exp(b_tm * tm[i]);
    qrs[4] = MDv[age] * exp(b_tm * tm[i]);
    qrs[5] = MWv[age] * exp(b_tm * tm[i]);
    qrs[6] = URv[age];
  };

  matrix<Type> Pm(N_Q, N_Q), Qm(N_Q, N_Q), cumPm(N_Q, N_Q);
  Type ll_val = 0;

  for (int i = 0; i < A.size(); i++)
  {
    if (fit[i] == 0)
    {
      fill_qrs(istart[i], i);
      Qm = KM(qrs, false, true);
      ll_val = log(Qm(A[i], Z[i]) + eps);
    }
    else if (fit[i] == 1)
    {
      fill_qrs(istart[i], i);
      Pm = KM(qrs, false, false);
      ll_val = log(Pm(A[i], Z[i]) + eps);
    }
    else if (fit[i] == 2)
    {
      cumPm = matrix<Type>::Identity(N_Q, N_Q);
      for (int t = istart[i]; t < end[i]; t++)
      {
        fill_qrs(t, i);
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
      int age = CppAD::Integer(sim_data(i, 0));
      qrs[0] = VXv[age];
      qrs[1] = VMv[age];
      qrs[2] = XMv[age] * exp(b_tx * sim_data(i, 1));
      qrs[3] = MSv[age] * exp(b_tm * sim_data(i, 2));
      qrs[4] = MDv[age] * exp(b_tm * sim_data(i, 2));
      qrs[5] = MWv[age] * exp(b_tm * sim_data(i, 2));
      qrs[6] = URv[age];
      PP.col(i) = KM(qrs, false, false);
    }
   
    REPORT(PP);
  }
  return dll;
}

