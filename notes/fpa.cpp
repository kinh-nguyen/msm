#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // Data
  DATA_VECTOR(y); // Response vector
  DATA_VECTOR(x); // Predictor vector

  // Parameters
  PARAMETER(beta0);      // Intercept
  PARAMETER(beta1);      // Coefficient for basis function 1
  PARAMETER(beta2);      // Coefficient for basis function 2
  PARAMETER(gamma);      // DFP-2 dynamic parameter 1
  PARAMETER(delta);      // DFP-2 dynamic parameter 2
  PARAMETER(log_sigma);  // Log of the residual standard deviation

  // Transformations
  Type sigma = exp(log_sigma);

  // Construct the 2x2 dynamics matrix A
  matrix<Type> A(2, 2);
  A(0, 0) = 0.0;
  A(0, 1) = 1.0;
  A(1, 0) = -delta;
  A(1, 1) = -gamma;

  // Initialize the negative log-likelihood
  Type nll = 0.0;

  // Loop through each observation
  for (int i = 0; i < x.size(); i++) {
    Type t_i = log(x(i));
    
    // Calculate the matrix exponential S = exp(A*t)
    // This is the core of the DFP-2 model
    matrix<Type> S_i = atomic::expm(A * t_i);

    // The basis functions are the top row of the solution matrix S
    Type B1_i = S_i(0, 0);
    Type B2_i = S_i(0, 1);

    // Calculate the linear predictor eta
    Type eta_i = beta0 + beta1 * B1_i + beta2 * B2_i;

    // Add the contribution to the negative log-likelihood
    nll -= dnorm(y(i), eta_i, sigma, true);
  }

  return nll;
}