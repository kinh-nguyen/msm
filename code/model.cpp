#include <TMB.hpp>

#define AGE_MAX 50
#define N_AGE 51
#define N_D 51 // differences, more than needed
#define N_PAR 7
#define N_Q 7
#define N_CC 37

using Eigen::seqN;

template <class T> 
struct Kube {
    vector<T> 
        masterQ = vector<T>(N_AGE * N_Q * N_Q * N_CC),
        masterP = vector<T>(N_AGE * N_Q * N_Q * N_CC),
        masterD = vector<T>(N_AGE * N_Q * N_Q * N_D * N_CC);
    matrix<T> 
        qM = matrix<T>(N_Q, N_Q),
        pM = matrix<T>(N_Q, N_Q), 
        pD = matrix<T>(N_Q, N_Q), 
        est = matrix<T>(N_AGE * N_CC, N_PAR);
    int len = N_Q * N_Q;
    Kube () {};
    Kube (matrix<T> modelmatrix, vector<T> betas)
    {
        masterP.setZero();
        masterQ.setZero();
        masterD.setZero();
        matrix<T> beta_i(N_CC + 2, 1);
        for (int i = 0; i < N_PAR; i++) {
            beta_i << betas({i, N_PAR + i}), betas(seqN(2*N_PAR + i*N_CC, N_CC));
            est.col(i) = modelmatrix * beta_i;
        }
        est = exp(est.array());
        for (int a = 0; a < N_AGE; a++) {
        for (int c = 0; c < N_CC; c++) {
            qM.setZero(); 
            int r = a*N_CC + c;
            qM(0, 1) = est(r, 0); // debut
            qM(0, 2) = est(r, 1); // marriage from virgin
            qM(1, 2) = est(r, 2); // marriage from debut
            qM(2, {3,4,5}) = est(r, {3,4,5}); // marriage dissolution
            qM({3,4,5}, 6) = est(r, {6,6,6}); // disso > remarried
            qM(6, {3,4,5}) = est(r, {3,4,5}); // remarried > disso = married > disso, we could add a(three) scaling parameter as well?
            qM.diagonal() = T(-1) * qM.rowwise().sum();
            memcpy(&masterQ(0) + a*N_CC*len + c*len, &qM(0), sizeof(T)*len);
            pM = expm(qM);
            memcpy(&masterP(0) + a*N_CC*len + c*len, &pM(0), sizeof(T)*len);
            for (int d = 0; d < N_D; d++) { // prep delta
                matrix<T> tmp = qM * d;
                pD = expm(tmp);
                memcpy(&masterD(0) + a*N_CC*N_D*len + c*N_D*len + d*len, &pD(0), sizeof(T)*len);
            }
            }
        }
    };
    matrix<T> operator()(int c, int a){
        matrix<T> ans = masterP(seqN(a*N_CC*len + c*len, len)).reshaped(N_Q, N_Q);
        return ans;
    };
    matrix<T> operator()(int c, int a, int d){
        matrix<T> ans = masterD(seqN(a*N_CC*N_D*len + c*N_D*len + d*len, len)).reshaped(N_Q, N_Q);
        return ans;
    }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);
  Type prior = 0.0;

  // data
  DATA_IVECTOR(afs);
  DATA_IVECTOR(aam);

  DATA_IVECTOR(VV);
  DATA_IVECTOR(VX);
  DATA_IVECTOR(VM);
  DATA_IVECTOR(XX);
  DATA_IVECTOR(XM);
  DATA_IVECTOR(MM);
  DATA_IVECTOR(MJ);
  DATA_IVECTOR(J);
  DATA_IVECTOR(delta);
  DATA_VECTOR(count);
  DATA_IVECTOR(cid);
  
  DATA_SPARSE_MATRIX(modelmatrix);

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
  for (int i = 0; i < N_PAR; i++) {
    vector<Type> f_cc = betas(seqN(N_PAR*2 + i*N_CC, N_CC));
    prior -= dnorm(f_cc, prior_cc(0), prior_cc(1), true).sum();
  }

  Kube<Type> KM(modelmatrix, betas);

  matrix<Type> o(afs.size(), 7);
  o.setOnes();

  for (int i = 0; i < afs.size(); i++) {
    if (VV[i] > 0)
        for (int j = 0; j <= VV[i]; j++)
            o(i, 0) *= KM(cid[i], j)(0,0);
    if (VX[i] > 0)
            o(i, 1) *= KM(cid[i], VX[i])(0,1);
    if (XX[i] > 0)
        for (int j = afs[i]; j <= afs[i] + XX[i]; j++) 
            o(i, 2) *= KM(cid[i], j)(1,1);
    if (XM[i] > 0)
            o(i, 3) *= KM(cid[i], XM[i])(1,2);
    if (MM[i] > 0)
        for (int j = aam[i]; j <= aam[i] + MM[i]; j++) 
            o(i, 4) *= KM(cid[i], j)(2,2);
    if (MJ[i] > 0)
            o(i, 5) *= KM(cid[i], MJ[i], delta[i])(2,J[i]);
    if (VM[i] > 0)
            o(i, 6) *= KM(cid[i], VM[i])(0,2);
  }
  REPORT(o);
  o = log(o.array()); // all the 1s disapear here
  vector<Type> oo = o.rowwise().sum();
  oo *= count; // element-wise
  dll += prior - oo.sum();

  REPORT(betas);
  REPORT(KM.est);
  REPORT(KM.masterP);
  REPORT(KM.masterQ);
  REPORT(KM.masterD);
  return dll;
}
