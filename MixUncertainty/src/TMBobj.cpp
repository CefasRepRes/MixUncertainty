#define TMB_LIB_INIT R_init_MixUncertainty
#include <TMB.hpp>
#include <contrib/OSA_multivariate_dists-main/distr.hpp>

// find cases of missing values
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// logit transformation
template<class Type>
vector<Type> transLogit(vector<Type> alpha){
  int dim = alpha.size();
  vector<Type> p(dim+1);
  vector<Type> expa = exp(alpha);
  Type s = sum(expa);
  Type lastp = Type(1);
  for(int i = 0; i < dim; ++i){
    p(i) = expa(i)/(Type(1) + s);
    lastp -= p(i);
  }
  p(dim)=lastp;
  return p;
}

// define objective function
template<class Type>
Type objective_function<Type>::operator()(){

  DATA_STRING(code);

  // Code definition:
  //   A = log-scale AR1 latent process with Multivariate Gaussian intervals,
  //       Gaussian observation error
  //   B = log-scale AR1 latent process with Gaussian intervals, Gaussian
  //       observation error
  //   C = multinomial logit-scale random walk latent process with Multivariate
  //       Gaussian intervals. Dirichlet observation error.
  //   D = multinomial logit-scale random walk latent process with Gaussian
  //       intervals. Dirichlet observation error.
  //   E = multinomial logit-scale AR1 latent process with Multivariate Gaussian
  //       intervals. Hurdle observation error using Dirichlet and Bernoulli.
  //   F = multinomial logit-scale AR1 latent process with Gaussian
  //       intervals. Hurdle observation error using Dirichlet and Bernoulli.

  if (code == "A") {
    #include "N_logMVN_rw.h"

  } else if (code == "B") {
    #include "N_logN_rw.h"

  } else if (code == "C") {
    #include "Dir_logitMVN_rw.h"

  } else if (code == "D") {
    #include "Dir_logitN_rw.h"

  } else if (code == "E") {
    #include "Dir_MVN_AR1_hurdle.h"

  } else if (code == "F") {
    #include "Dir_N_AR1_hurdle.h"

  } else {
    error("Unknown model type")
  }
  return 0;
}
