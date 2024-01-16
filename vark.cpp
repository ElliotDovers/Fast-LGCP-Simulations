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
   PARAMETER_VECTOR(mu);         // Random effect q-dist means - length K
   PARAMETER_VECTOR(logsig2);    // Random effect q-dist sd - length K

   ////////////////////////////////////////////////////////////////////
   ///// Calculate the likelihood at the presence and quad points /////
   ////////////////////////////////////////////////////////////////////  
   
   // pre-calc some matrix computations - uses sparse matrices
   vector<Type> Beta(2);
   Beta[0] = b0;
   Beta[1] = b1;
   vector<Type> sig2 = exp(logsig2);
   vector<Type> XB_pres = X_pres * Beta;
   vector<Type> XB_quad = X_quad * Beta;
   vector<Type> Zmu_pres = Z_pres * mu;
   vector<Type> Zmu_quad = Z_quad * mu;
   SparseMatrix<Type> Z2= Z_quad.cwiseProduct(Z_quad);
   vector<Type> ZsigZ = Z2 * sig2;
      
   Type ll1 = 0.0; // initialise the log-likelihood component 1
   Type ll2 = 0.0; // initialise the log-likelihood component 2
   ll1 += XB_pres.sum() + Zmu_pres.sum();
   ll2 -= (quad_size * exp(XB_quad + Zmu_quad + (0.5 * ZsigZ))).sum();
   
   /////////////////////////////////////////////////////////////////////////
   //////////// KL Divergence between q() and p() - Component 3 ////////////
   /////////////////////////////////////////////////////////////////////////
   
   // Obtain MLE estimate of the prior variance - via profiling
   Type del2;
   vector<Type> tmp(Z_pres.cols());
   for (int k = 0; k < Z_pres.cols(); ++k) {
      tmp(k) = sig2(k) + (mu(k) * mu(k));
   }
   del2 = tmp.sum() / tmp.size();
   Type DKL = (Z_pres.cols() * log(del2)) - logsig2.sum();
   Type ll3 = -0.5 * DKL;
   
   /////////////////////////////////////////////////////////////////////////
   
   // Combine for the loglikelihood
   Type nll = -(ll1 + ll2 + ll3);
   // ADREPORT(Beta);
   ADREPORT(sig2);
   ADREPORT(del2);
   REPORT(ll1);
   REPORT(ll2);
   REPORT(ll3);
      
   return nll;
}
