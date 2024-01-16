// log-likelihood for presence-only data Inhomogeneous Poisson Process - Fixed Rank Kriging Latent Approximation
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
      DATA_MATRIX(X_pres);          // Fixed effect (N+M)xP matrix - at the pres & quad pts
      DATA_MATRIX(X_quad);          // Fixed effect (N+M)xP matrix - at the pres & quad pts
      DATA_VECTOR(quad_size);       // Vector describing the quadrature sizes - length N + M
      PARAMETER(b0);             // Fixed effect Intercept
      PARAMETER(b1);             // Fixed effect Slope
      
      Type ll1 = 0.0; // initialise the log-likelihood component 1
      Type ll2 = 0.0; // initialise the log-likelihood component 2
      
      ////////////////////////////////////////////////////////////////
      // Conditional log-likelihood at the presence and quad points //
      ////////////////////////////////////////////////////////////////
      
      // Initial elements to required for later calculations
      vector<Type> Beta(2);
      Beta[0] = b0;
      Beta[1] = b1;
      vector<Type> XB_pres = X_pres * Beta;
      vector<Type> XB_quad = X_quad * Beta;
      
      ll1 += XB_pres.sum();
      ll2 -= (quad_size * exp(XB_quad)).sum();
      
      // Combine for the loglikelihood
      Type nll = -(ll1 + ll2);
      REPORT(ll1);
      REPORT(ll2);
      
      return nll;
      }
