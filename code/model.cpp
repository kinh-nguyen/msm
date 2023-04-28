#include <TMB.hpp>

#define AGE_MAX 50
#define N_AGE 51
#define N_D 51 // differences, more than needed
#define N_PAR 8
#define N_Q 7
#define N_YOB 37 // 1965-2001

using Eigen::seqN;

template <class T> 
struct Kube {
    vector<T>
        intercept = vector<T>(N_PAR),
        beta_t = vector<T>(N_PAR),
        yobsm = vector<T>(N_YOB * N_PAR),
        q_v = vector<T>(N_PAR);
    matrix<T>
        qM = matrix<T>(N_Q, N_Q),
        pM = matrix<T>(N_Q, N_Q);
    Kube () {};
    Kube (vector<T> betas, vector<T> yobsm) :
        intercept(betas(seqN(0, N_PAR))),
        beta_t   (betas(seqN(N_PAR, N_PAR))),
        yobsm    (yobsm)
    {};
    matrix<T> operator()(int a, int y, T delta = 1) {
        for (int i = 0; i < N_PAR; i++)
            q_v(i) = intercept(i) + beta_t(i) * a + yobsm(i * N_YOB + y);
        q_v = exp(q_v);
        qM.setZero(); 
        qM(0, 1)       = q_v(0); // debut
        qM(1, 2)       = q_v(1); // marriage
        qM(2, {3,4,5}) = q_v({2,3,4}); // marriage dissolution
        qM({3,4,5}, 6) = q_v({5,6,7}); // disso > remarried
        qM(6, {3,4,5}) = q_v({2,3,4}); // remarried > disso = married > disso, we could add a(three) scaling parameter as well?
        qM.diagonal() = T(-1) * qM.rowwise().sum();
        matrix<T> tmp = qM * delta;
        pM = expm(qM);
        return pM;
    };
};

template <class T> // 3 or 4 need to implement tophi
struct mAR {
    matrix<T> Sigma;
    vector<T> phi;
    mAR () {};
    mAR (int m) : Sigma(m, m), phi(m) {Sigma.setIdentity();};
    vector<T> to_phi (vector<T> x) {
        if (x.size() != 2) Rf_error("not implemented");
        vector<T> psi(2);
        psi[0] = 2. * exp(x[0]) / (1. + exp(x[0])) - 1.;
        psi[1] = 2. * exp(x[1]) / (1. + exp(x[1])) - 1.;
        phi[1] = psi[1];
        phi[0] = psi[0] * (1.0 - phi[1]);
        if (phi[1] == -1) phi[1] += FLT_EPSILON;
        if (phi[1] == 1 - phi.abs()(0) ) phi[1] -= DBL_EPSILON;
        return phi;
    };
    T operator()(vector<T> pacf, vector<T> x, bool zeroing = true){
        T dll = 0.;
        dll += density::MVNORM(Sigma)(pacf);
        phi = to_phi(pacf);
        dll += density::ARk(phi)(x);
        if (zeroing)
            dll -= dnorm(sum(x), T(0), T(0.001) * x.size(), true);
        return dll;
    };
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);
  Type prior = 0.0;

  // data
  DATA_IVECTOR(afs);
  DATA_IVECTOR(aam);
  DATA_IVECTOR(yob);

  DATA_IVECTOR(VV);
  DATA_IVECTOR(VX);
  DATA_IVECTOR(XX);
  DATA_IVECTOR(XM);
  DATA_IVECTOR(MM);
  DATA_IVECTOR(MJ);
  DATA_IVECTOR(J);
  DATA_IVECTOR(delta);
  DATA_VECTOR(w);

  // priors
  DATA_VECTOR(sd_q);

  // Coefs
  PARAMETER_VECTOR(betas); 
  // base rate | intercept
  vector<Type> intercepts = betas(seqN(0, N_PAR));
  prior -= dnorm(intercepts, sd_q(0), sd_q(1), true).sum();
  // Soft-constraints remarried to be the same
  prior -= dnorm(betas(5) - betas(6), Type(0), Type(0.001), true);
  prior -= dnorm(betas(5) - betas(7), Type(0), Type(0.001), true);
  // age's coeff | time in the hazard
  vector<Type> beta_t = betas(seqN(N_PAR, N_PAR));
  prior -= dnorm(beta_t, Type(0), Type(0.5), true).sum();
  // Soft-constraints remarried to be the same
  prior -= dnorm(betas(N_PAR + 5) - betas(N_PAR + 6), Type(0), Type(0.001), true);
  prior -= dnorm(betas(N_PAR + 5) - betas(N_PAR + 7), Type(0), Type(0.001), true);
  // Temporal on year of births
  PARAMETER_VECTOR(yobsm);
  PARAMETER_VECTOR(pacf);
  //AR2 prior for each term
  mAR<Type> MAR2(2);
  for (int i = 0; i < N_PAR; i++)
    prior += MAR2(pacf, yobsm(seqN(i * N_YOB, N_YOB)));
  // Soft-constraints effects yob on remarried to be the same
  vector<Type> 
  diff = yobsm(seqN(5*N_YOB, N_YOB)) - yobsm(seqN(6*N_YOB, N_YOB));
  prior -= dnorm(diff, Type(0), Type(0.001 * N_YOB), true).sum();
  diff = yobsm(seqN(5*N_YOB, N_YOB)) - yobsm(seqN(7*N_YOB, N_YOB));
  prior -= dnorm(diff, Type(0), Type(0.001 * N_YOB), true).sum();

  Kube<Type> KM(betas, yobsm);

  vector<Type> o(afs.size());
  o.setOnes();
  for (int i = 0; i < afs.size(); i++) {
    if (VV[i] > 0)
        for (int j = 0; j <= VV[i]; j++)
            o(i) *= KM(j, yob(i))(0,0);
    if (VX[i] > 0)
            o(i) *= KM(VX[i], yob(i))(0,1);
    if (XX[i] > 0)
        for (int j = afs[i]; j <= afs[i] + XX[i]; j++) 
            o(i) *= KM(j, yob(i))(1,1);
    if (XM[i] > 0)
            o(i) *= KM(XM[i], yob(i))(1,2);
    if (MM[i] > 0)
        for (int j = aam[i]; j <= aam[i] + MM[i]; j++) 
            o(i) *= KM(j, yob(i))(2,2);
    if (MJ[i] > 0)
            o(i) *= KM(MJ[i], yob(i), delta[i])(2,J[i]);
  }
  REPORT(o);
  vector<Type> tmp = log(o) * w.array(); 
  dll += prior - tmp.sum();

  SIMULATE {
    matrix<Type> tmp2(N_Q, N_Q);
    int len = N_Q * N_Q;
    vector<Type> estP(len * N_AGE * N_YOB);
    for (int j = 0; j < N_YOB; j++)
      for (int i = 0; i < N_AGE; i++) {
        tmp2 = KM(i, j);
        memcpy(&estP(0) + i*len + j*N_AGE*len, &tmp2(0), sizeof(Type) * len);
      }
    REPORT(estP);
  }
  REPORT(tmp);
  REPORT(betas);
  REPORT(yobsm);
  REPORT(pacf);
  return dll;
}
