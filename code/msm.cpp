#include <TMB.hpp>
#include "ktools.hpp"

#define N_PAR 7
#define N_Q 7

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

// Log-logistic baseline hazard function where
// t is time, exp_b0 = alpha, b1 is beta > 0
template <class Type>
Type LLGx(Type t, Type exp_b0, Type b1) {
  if (t == 0.0) return Type(0.0); 
  Type 
    t_pow_b1 = pow(t, b1),
    numerator = exp_b0 * b1 * t_pow_b1 / t,
    denominator = Type(1.0) + exp_b0 * t_pow_b1;
  return numerator / denominator;
}

// Log-skew-logistic baseline hazard function 
template <class Type>
Type LSLx(Type t, Type lamda, Type p, Type gamma) {
  if (t == 0.0) return Type(0.0); // can be undefined without reparameterization
  Type
      u = pow(lamda * t, -p),
      n = gamma * p * u,
      d = t * (1 + u) * (pow(1 + u, gamma) - Type(1.0));
  return n / d;
}

template<class Type>
Type EDLLx(
    Type t, Type log_t,         // Current time t and log(t)
    Type log_alpha,             // log(alpha) - LL scale related (alpha > 0)
    Type beta,                  // beta - LL shape (beta > 0)
    Type log_beta,              // log(beta)
    Type T_M,                   // T_M - modifier timing (can be any real number)
    Type k                      // k - modifier rate/steepness (k > 0)
) {
    Type log_t_div_alpha = log_t - log_alpha;
    // log(h_LL(t)) = log(beta) - log(alpha) + (beta-1)*log(t/alpha) - log(1 + (t/alpha)^beta)
    Type log_LL_denom_factor = logspace_add(Type(0.0), beta * log_t_div_alpha);
    Type log_h_LL = log_beta - log_alpha + (beta - Type(1.0)) * log_t_div_alpha - log_LL_denom_factor;
    // Modifier part: log(1 + exp(k*(t-T_M)))
    Type log_modifier_denom = logspace_add(Type(0.0), k * (t - T_M));
    Type log_final_hz = log_h_LL - log_modifier_denom;
    return exp(log_final_hz);
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
  vector<Type> log_t = log(dstart);
  DATA_IVECTOR(end);
  DATA_VECTOR(n);
  DATA_IVECTOR(fit);
  DATA_INTEGER(n_age);

  PARAMETER_VECTOR(log_T_M);
  vector<Type> T_M = exp(log_T_M);
  dll -= dnorm(log_T_M, Type(0), Type(5), true).sum();

  PARAMETER_VECTOR(log_beta);
  vector<Type> beta = exp(log_beta);
  dll -= dnorm(log_beta, Type(0), Type(5), true).sum();

  PARAMETER_VECTOR(log_alpha);
  vector<Type> alpha = exp(log_alpha);
  dll -= dnorm(log_alpha, Type(0), Type(5), true).sum();

  PARAMETER_VECTOR(log_k);
  vector<Type> k = exp(log_k);
  dll -= dnorm(log_k, Type(0), Type(5), true).sum();

  PQ<Type> KM;
  vector<Type> q_rs(N_PAR);

  vector<Type> indiv_ll(A.size());
  indiv_ll.setZero();

  auto fill_qrs = [&](Type age, Type log_age) {
    for (int p = 0; p < N_PAR; p++)
      q_rs[p] = EDLLx(age, log_age, log_alpha[p], beta[p], log_beta[p], T_M[p], k[p]);
  };

  for (int i = 0; i < A.size(); i++) {
    Type ll = 0;
    fill_qrs(dstart[i], log_t[i]);
    if (fit[i] == 0) 
      ll = -n[i] * KM(q_rs, true, true)(A[i], Z[i]); 
    else if (fit[i] == 1) 
      ll = -n[i] * KM(q_rs, true, false)(A[i], Z[i]); 
    else if (fit[i] == 2) {
      matrix<Type> cumP = matrix<Type>::Identity(N_Q, N_Q);
      for (int t = istart[i]; t < end[i]; t++)
      {
        fill_qrs(Type(t), log(Type(t)));
        cumP = cumP * KM(q_rs, false, false);
      }
      ll = -n[i] * log(cumP(A[i], Z[i]));
    }
    dll += ll;
    indiv_ll[i] = ll;
  }
  
  SIMULATE {
    array<Type> qM_rep(N_Q, N_Q, n_age), pM_rep(N_Q, N_Q, n_age);
    for (int s = 1; s < n_age; ++s) {
      fill_qrs(Type(s), log(Type(s)));
      matrix<Type> qM_mat = KM(q_rs, false, true);
      matrix<Type> pM_mat = KM(q_rs, false, false);
      for (int i = 0; i < N_Q; ++i)
        for (int j = 0; j < N_Q; ++j) {
          qM_rep(i, j, s) = qM_mat(i, j);
          pM_rep(i, j, s) = pM_mat(i, j);
        }
    }
    REPORT(qM_rep);
    REPORT(pM_rep);
    REPORT(T_M);
    REPORT(beta);
    REPORT(alpha);
    REPORT(k);
    REPORT(indiv_ll);
  }
  return dll;
}
