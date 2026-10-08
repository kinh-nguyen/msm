#include <TMB.hpp>

#define c_const 0.8

template <class Type>
struct gk_quantile
{
  Type A, B, g, k, t;
  // Type c_const = 0.8;

  gk_quantile(Type A, Type B, Type g, Type k, Type t)
      : A(A), B(B), g(g), k(k), t(t) {}

  Type operator()(vector<Type> z)
  {
    Type z0 = z[0];
    if (abs(z0) > 20.0)
      return Type(1e10);
    Type tanh_term = tanh(g * z0 / Type(2.0));
    Type base = 1.0 + z0 * z0;
    Type base_k = pow(base, k);
    // Type base = pow(Type(1.0) + z0 * z0, k);
    // Type Q = A + B * (Type(1.0) + c_const * tanh_term) * base * z0;
    Type Q = A + B * (1.0 + c_const * tanh_term) * base_k * z0;
    return Q - t;
  }
};

template <class Type>
Type objective_function<Type>::operator()()
{
  using newton::Newton;
  using newton::newton_config_t;

  // Data
  DATA_SCALAR(t);                    // observed value
  DATA_STRUCT(cfg, newton_config_t); // Newton config

  // Parameters
  PARAMETER(A);
  PARAMETER(B);
  PARAMETER(g);
  PARAMETER(k);

  // Solve Q(z) = t
  gk_quantile<TMBad::ad_aug> F(A, B, g, k, t);
  vector<Type> start(1);
  start << 0.0;
  vector<Type> z = Newton(F, start, cfg);

  // Derivative Q'(z)
  Type z0 = z[0];
  Type tanh_term = tanh(g * z0 / Type(2.0));
  Type sech2 = 1.0 - tanh_term * tanh_term; // More stable than 1/cosh²(x)
  Type base = 1.0 + z0 * z0;

  // For efficiency, call pow() once and derive the second power term.
  Type base_k = pow(base, k);
  Type base_k_minus_1 = base_k / base;

  // The derivative using the product rule: B * [f(z)h'(z) + f'(z)h(z)]
  Type dQdz =
      B * (1.0 + c_const * tanh_term) * (base_k + 2.0 * k * z0 * z0 * base_k_minus_1) +
      B * c_const * (g / 2.0) * sech2 * base_k * z0;

  // Compute hazard function
  Type phi = dnorm(z0, Type(0.0), Type(1.0), false);
  Type S = 1.0 - pnorm(z0);
  Type hazard = phi / (S * dQdz);

  REPORT(z);
  REPORT(hazard);

  return -log(hazard); // Example: negative log-likelihood
}
