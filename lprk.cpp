// Variational Approx. of log-likelihood for presence-only data Cox Process
// Single Resolution, non-spatially correlated random coefficients
// Sparse Basis Functions 
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
   
   using namespace Eigen;
   
   DATA_MATRIX(X_pres);          // Fixed effect (N+M)xP matrix - at the pres & quad pts
   DATA_SPARSE_MATRIX(Z_pres);   // Basis functions (N+M)xK matrix - at the pres & quad pts
   DATA_MATRIX(X_quad);          // Fixed effect (N+M)xP matrix - at the pres & quad pts
   DATA_SPARSE_MATRIX(Z_quad);   // Basis functions (N+M)xK matrix - at the pres & quad pts
   DATA_VECTOR(quad_size);       // Vector describing the quadrature sizes - length N + M
   PARAMETER(b0);             // Fixed effect Intercept
   PARAMETER(b1);             // Fixed effect Slope
   PARAMETER_VECTOR(u);          // Random effect q-dist means - length K
   PARAMETER(logdel);           // Random coefficient variance

   vector<Type> Beta(2);
   Beta[0] = b0;
   Beta[1] = b1;
   Type ll1 = 0.0; // initialise the log-likelihood component 1
   Type ll2 = 0.0; // initialise the log-likelihood component 2
   Type ll3 = 0.0; // initialise the log-likelihood component 3
   Type mu = 0.0; // initialise the random coefficient mean
   Type del = exp(logdel);
   
   // The random coefficient contribution
   for (int k = 0; k < Z_pres.cols(); ++k) {
      ll3 += dnorm(u(k), mu, del, true);
   }
   
   ////////////////////////////////////////////////////////////////////
   ///// Calculate the likelihood at the presence and quad points /////
   ////////////////////////////////////////////////////////////////////  
   
   // pre-calc some matrix computations - uses sparse matrices
   vector<Type> XB_pres = X_pres * Beta;
   vector<Type> XB_quad = X_quad * Beta;
   vector<Type> Zu_pres = Z_pres * u;
   vector<Type> Zu_quad = Z_quad * u;
   
   ll1 += XB_pres.sum() + Zu_pres.sum();
   ll2 -= (quad_size * exp(XB_quad + Zu_quad)).sum();
   
   /////////////////////////////////////////////////////////////////////////
   
   // Combine for the loglikelihood
   Type nll = -(ll1 + ll2 + ll3);
   Type del2 = del * del;
   // ADREPORT(Beta);
   ADREPORT(del2);
   REPORT(ll1);
   REPORT(ll2);
   REPORT(ll3);
      
   return nll;
}
