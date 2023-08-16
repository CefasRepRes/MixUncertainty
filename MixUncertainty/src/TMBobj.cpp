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

  // Read data and model type
  // ========================

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

  DATA_STRING(code);

  // read in data
  DATA_ARRAY(dat); // matrix: variable x years

  // Read parameters
  // ===============
  //
  // Latent process
  // --------------

  // Common parameters
  PARAMETER_ARRAY(rw);        // matrix of unobserved "true" values
  PARAMETER_VECTOR(logSdMVN); // variance - MVN

  // Model-specific parameters
  PARAMETER(logitRho);           // Multivariate Normal latent correlation parameter
  PARAMETER_VECTOR(mu);          // stationarity mean for each element
  PARAMETER_VECTOR(logitARrho);  // rho parameter (correlation) for autoregressive process (AR1)

  // Observation process
  // -------------------

  // Univariate Normal
  PARAMETER_VECTOR(logSdObs);    // Univariate Normal observation variance

  // Dirichlet
  PARAMETER(logTau); // concentration parameter for Dirichlet observation

  // Bernoulli
  PARAMETER(logb0); // linear predictor intercept
  PARAMETER(logb1); // linear predictor coefficient

  // Missing observations
  PARAMETER_VECTOR(dat_missing); // vector of missing values

  // Transform parameters
  // --------------------

  // data dimensions
  int nyear = rw.dim(1);
  int nvar  = rw.dim(0);

  // transform parameters
  vector <Type> SdMVN = exp(logSdMVN);
  vector <Type> SdObs = exp(logSdObs);
  Type rho = Type(2) * exp(logitRho)/(exp(logitRho) + Type(1)) - Type(1);
  vector <Type> ARrho = Type(2) * exp(logitARrho)/(exp(logitARrho) + Type(1)) - Type(1);
  vector <Type> sqRho = sqrt(Type(1) - pow(ARrho, Type(2)));
  Type            tau = exp(logTau);
  Type             b0 = exp(logb0) * Type(-1);
  Type             b1 = exp(logb1);

  // Implement appropriate model
  // ===========================

  Type nll = 0;

  if (code == "A") {

    // ========================================================================================
    // Multivariate normal random walk on log-scale. Univariate normal observation on log-scale
    // ========================================================================================

    // loop over each matrix element and identify missing parameters
    int missing_count = 0;
    for (int i = 0; i < nvar; i++) {
      for (int j = 0; j < nyear; j++) {
        if(isNA(dat(i,j))) {
          dat(i,j) = dat_missing(missing_count);
          missing_count++;
        }
      }
    }

    // prepare variance-covariance matrix
    // ----------------------------------

    matrix <Type> MVNcov(nvar, nvar);
    MVNcov.setZero();
    for(int i = 0; i < nvar; i++) {
      for(int j = 0; j < nvar; j++) {
        if(i == j) {
          MVNcov(i,j) = pow(SdMVN(i), 2); // variances
        } else {
          MVNcov(i,j) = SdMVN(i) * SdMVN(j) * rho; // covariances
        }
      }
    }

    // Compute negative log-likelihood
    // -------------------------------

    using namespace density;
    MVNORM_t<Type> neg_log_density(MVNcov);

    // calculate density at start of AR1 process
    nll += neg_log_density(rw.col(0) - mu);

    // Stationary AR process on the expectation
    for (int i = 1; i < nyear; ++i){
      nll += neg_log_density((rw.col(i) - (mu + ARrho * (rw.col(i-1) - mu)))/sqRho);
    }

    // Observation process
    for (int i = 0; i < nyear; ++i){
      nll -= dnorm(vector<Type>(dat.col(i)),
                   vector<Type>(rw.col(i)),
                   SdObs, true).sum();
    }

    // report fitted variance-covariance matrix
    REPORT(MVNcov);
    REPORT(rho);
    REPORT(dat_missing);

    // report random parameters with standard error
    REPORT(rw);

  } else if (code == "B") {

    // ======================================================================================
    // Univariate normal random walk on log-scale. Univariate normal observation on log-scale
    // ======================================================================================

    // loop over each matrix element and identify missing parameters
    int missing_count = 0;
    for (int i = 0; i < nyear; i++) {
      if(isNA(dat(0,i))) {
        dat(0,i) = dat_missing(missing_count);
        missing_count++;
      }
    }

    // prepare variance-covariance matrix
    // ----------------------------------

    matrix <Type> MVNcov(nvar, nvar);
    MVNcov.setZero();
    for(int i = 0; i < nvar; i++) {
          MVNcov(i,i) = pow(SdMVN(i), 2); // variances
    }

    // Compute negative log-likelihood
    // -------------------------------

    using namespace density;
    MVNORM_t<Type> neg_log_density(MVNcov);

    // calculate density at start of AR1 process
    nll += neg_log_density(rw.col(0) - mu);

    // Stationary AR process on the expectation
    for (int i = 1; i < nyear; ++i){
      nll += neg_log_density((rw.col(i) - (mu + ARrho * (rw.col(i-1) - mu)))/sqRho);
    }

    // Observation process
    for (int i = 0; i < nyear; ++i){
      Type rwi  = rw(0,i);
      Type dati = dat(0,i);

      nll -= dnorm(dati, rwi, SdObs(0), true);
    }

    // report variance
    REPORT(MVNcov);

    // report fitted missing values
    REPORT(dat_missing);

    // report random parameters with standard error
    REPORT(rw);

  } else if (code == "C") {

    // =============================================================================
    // Multivariate normal random walk on logit-scale. Dirichlet observation
    // =============================================================================

    // prepare matrix to store alpha values
    // ------------------------------------

    matrix <Type> alphaMatrix(nvar+1, nyear);
    alphaMatrix.setZero();

    // prepare variance-covariance matrix
    // ----------------------------------

    matrix <Type> MVNcov(nvar, nvar);
    MVNcov.setZero();
    for(int i = 0; i < nvar; i++) {
      for(int j = 0; j < nvar; j++) {
        if(i == j) {
          MVNcov(i,j) = pow(SdMVN(i), 2); // variances
        } else {
          MVNcov(i,j) = SdMVN(i) * SdMVN(j) * rho; // covariances
        }
      }
    }

    // Compute negative log-likelihood
    // -------------------------------

    using namespace density;
    MVNORM_t<Type> neg_log_density(MVNcov);

    // Random walk process on the expectation
    for (int i = 1; i < nyear; ++i){
      nll += neg_log_density(rw.col(i) - rw.col(i-1));
    }

    // Dirichlet observation process
    for (int i = 0; i < nyear; ++i){

      // extract i'th compositional data vector
      vector <Type> datVector = dat.col(i);

      // calculate i'th alpha vector
      vector <Type> betaVector = rw.col(i);
      vector <Type> alphaVector = tau * transLogit(betaVector);

      // store alpha vector
      alphaMatrix.col(i) = alphaVector;

      // find number of non-zero values
      int true_count = 0;
      for (int j = 0; j < datVector.size(); ++j) {
        if (datVector(j) > Type(0)) {
          true_count++;
        }
      }

      // only compute observation likelihood if data vector has more than 1 non-zero element
      if (true_count > 1) {

        // generate index of non-zero values
        vector <int> datIdx(true_count);
        true_count = 0;
        for (int j = 0; j < datVector.size(); ++j) {
          if (datVector(j) > Type(0)) {
            datIdx(true_count) = j;
            true_count++;
          }
        }

        // extract non-zero elements
        vector <Type> datNonZero = datVector(datIdx);
        vector <Type> alphaNonZero = alphaVector(datIdx);

        // negative log-likelihood
        nll += -ddirichlet(datNonZero,
                           alphaNonZero,
                           true);
      }
    }

    // report fitted variance-covariance matrix
    REPORT(MVNcov);

    // report random parameters
    REPORT(alphaMatrix);

  } else if (code == "D") {

    // =============================================================================
    // Multivariate normal random walk on log-scale. Dirichlet observation
    // =============================================================================

    // prepare matrix to store alpha values
    // ------------------------------------

    matrix <Type> alphaMatrix(nvar+1, nyear);
    alphaMatrix.setZero();

    // prepare variance-covariance matrix
    // ----------------------------------

    matrix <Type> MVNcov(nvar, nvar);
    MVNcov.setZero();
    for(int i = 0; i < nvar; i++) {
      MVNcov(i,i) = pow(SdMVN(i), 2); // variances
    }

    // Compute negative log-likelihood
    // -------------------------------

    using namespace density;
    MVNORM_t<Type> neg_log_density(MVNcov);

    // Random walk process on the expectation
    for (int i = 1; i < nyear; ++i){
      nll += neg_log_density(rw.col(i) - rw.col(i-1));
    }

    // Dirichlet observation process
    for (int i = 0; i < nyear; ++i){

      // extract i'th compositional data vector
      vector <Type> datVector = dat.col(i);

      // calculate i'th alpha vector
      vector <Type> betaVector = rw.col(i);
      vector <Type> alphaVector = tau * transLogit(betaVector);

      // store alpha vector
      alphaMatrix.col(i) = alphaVector;

      // find number of non-zero values
      int true_count = 0;
      for (int j = 0; j < datVector.size(); ++j) {
        if (datVector(j) > Type(0)) {
          true_count++;
        }
      }

      // only compute observation likelihood if data vector has more than 1 non-zero element
      if (true_count > 1) {

        // generate index of non-zero values
        vector <int> datIdx(true_count);
        true_count = 0;
        for (int j = 0; j < datVector.size(); ++j) {
          if (datVector(j) > Type(0)) {
            datIdx(true_count) = j;
            true_count++;
          }
        }

        // extract non-zero elements
        vector <Type> datNonZero = datVector(datIdx);
        vector <Type> alphaNonZero = alphaVector(datIdx);

        // negative log-likelihood
        nll += -ddirichlet(datNonZero,
                           alphaNonZero,
                           true);
      }
    }

    // report fitted variance-covariance matrix
    REPORT(MVNcov);

    // report random parameters
    REPORT(alphaMatrix);

  } else if (code == "E") {

    // =============================================================================
    // Negative log-likelihood for multivariate normal AR1 on logit-scale (Hurdle)
    // =============================================================================

    // prepare matrix to store alpha values
    // ------------------------------------

    matrix <Type> alphaMatrix(nvar+1, nyear);
    alphaMatrix.setZero();

    // prepare matrix to store observation probability
    // -----------------------------------------------

    matrix <Type> probabilityMatrix(nvar+1, nyear);
    probabilityMatrix.setZero();

    // prepare variance-covariance matrix
    // ----------------------------------

    matrix <Type> MVNcov(nvar, nvar);
    MVNcov.setZero();
    for(int i = 0; i < nvar; i++) {
      for(int j = 0; j < nvar; j++) {
        if(i == j) {
          MVNcov(i,j) = pow(SdMVN(i), Type(2)); // variances
        } else {
          MVNcov(i,j) = SdMVN(i) * SdMVN(j) * rho; // covariances
        }
      }
    }

    // Compute negative log-likelihood
    // -------------------------------

    using namespace density;
    MVNORM_t<Type> neg_log_density(MVNcov);

    // calculate density at start of AR1 process
    nll += neg_log_density(rw.col(0) - mu);

    // Stationary AR process on the expectation
    for (int i = 1; i < nyear; ++i){
      nll += neg_log_density((rw.col(i) - (mu + ARrho * (rw.col(i-1) - mu)))/sqRho);
    }

    // Observation process
    for (int i = 0; i < nyear; ++i){

      // Calculate expectation for each compositional element
      // ----------------------------------------------------

      // extract i'th compositional data vector
      vector <Type> datVector = dat.col(i);

      // calculate i'th expectation vector
      vector <Type> rwVector          = rw.col(i);
      vector <Type> expectationVector = transLogit(rwVector);

      // Bernoulli likelihood
      // --------------------

      // how many zero elements in the vector?
      vector <Type> x(nvar+1);

      // define a vector with elements 1 if data <= 0 and 0 if data > 0
      for (int j = 0; j < nvar+1; j++){
        if (!isNA(datVector(j))) {
          if (datVector(j) <= 0.000000001) {
            x(j) = Type(0);
          } else {
            x(j) = Type(1);
          }

          // calculate linear predictor
          Type l = b0 + b1 * expectationVector(j);

          // logit transform linear predictor to p
          Type p = Type(1)/(Type(1) + exp(-l));

          // store probability value in matrix
          probabilityMatrix(j,i) = p;

          // calculate log-likelihood
          nll -= x(j) * log(p) + (1 - x(j)) * log(1 - p);
        }
      }

      // Dirichlet likelihood
      // --------------------

      // calculate alpha vector
      vector <Type> alphaVector = tau * expectationVector;

      // store alpha vector
      alphaMatrix.col(i) = alphaVector;

      // find number of non-zero values
      int true_count = 0;
      for (int j = 0; j < datVector.size(); ++j) {
        if (!isNA(datVector(j))) {
          if (datVector(j) > Type(0)) {
            true_count++;
          }
        }
      }

      // only compute observation likelihood if data vector has more than 1 non-zero element
      if (true_count > 1) {

        // generate index of non-zero values
        vector <int> datIdx(true_count);
        true_count = 0;
        for (int j = 0; j < datVector.size(); ++j) {
          if (!isNA(datVector(j))) {
            if (datVector(j) > Type(0)) {
              datIdx(true_count) = j;
              true_count++;
            }
          }
        }

        // extract non-zero elements
        vector <Type> datNonZero   = datVector(datIdx);
        vector <Type> alphaNonZero = alphaVector(datIdx);

        // negative log-likelihood
        nll += -ddirichlet(datNonZero,
                           alphaNonZero,
                           true);
      }
    }

    // report fitted variance-covariance matrix
    REPORT(MVNcov);

    // report random parameters
    REPORT(alphaMatrix);
    REPORT(probabilityMatrix);

  } else if (code == "F") {

    // =============================================================================
    // Negative log-likelihood for univariate normal AR1 on logit-scale (Hurdle)
    // =============================================================================

    // prepare matrix to store alpha values
    // ------------------------------------

    matrix <Type> alphaMatrix(nvar+1, nyear);
    alphaMatrix.setZero();

    // prepare matrix to store observation probability
    // -----------------------------------------------

    matrix <Type> probabilityMatrix(nvar+1, nyear);
    probabilityMatrix.setZero();

    // prepare variance-covariance matrix
    // ----------------------------------

    matrix <Type> MVNcov(nvar, nvar);
    MVNcov.setZero();
    for(int i = 0; i < nvar; i++) {
      MVNcov(i,i) = pow(SdMVN(i), Type(2)); // variances

    }

    // Compute negative log-likelihood
    // -------------------------------

    using namespace density;
    MVNORM_t<Type> neg_log_density(MVNcov);

    // calculate density at start of AR1 process
    nll += neg_log_density(rw.col(0) - mu);

    // Stationary AR process on the expectation
    for (int i = 1; i < nyear; ++i){
      nll += neg_log_density((rw.col(i) - (mu + ARrho * (rw.col(i-1) - mu)))/sqRho);
    }

    // Observation process
    for (int i = 0; i < nyear; ++i){

      // Calculate expectation for each compositional element
      // ----------------------------------------------------

      // extract i'th compositional data vector
      vector <Type> datVector = dat.col(i);

      // calculate i'th expectation vector
      vector <Type> rwVector          = rw.col(i);
      vector <Type> expectationVector = transLogit(rwVector);

      // Bernoulli likelihood
      // --------------------

      // how many zero elements in the vector?
      vector <Type> x(nvar+1);

      // define a vector with elements 1 if data <= 0 and 0 if data > 0
      for (int j = 0; j < nvar+1; j++){
        if (!isNA(datVector(j))) {
          if (datVector(j) <= 0.000000001) {
            x(j) = Type(0);
          } else {
            x(j) = Type(1);
          }

          // calculate linear predictor
          Type l = b0 + b1 * expectationVector(j);

          // logit transform linear predictor to p
          Type p = Type(1)/(Type(1) + exp(-l));

          // store probability value in matrix
          probabilityMatrix(j,i) = p;

          // calculate log-likelihood
          nll -= x(j) * log(p) + (1 - x(j)) * log(1 - p);
        }
      }

      // Dirichlet likelihood
      // --------------------

      // calculate alpha vector
      vector <Type> alphaVector = tau * expectationVector;

      // store alpha vector
      alphaMatrix.col(i) = alphaVector;

      // find number of non-zero values
      int true_count = 0;
      for (int j = 0; j < datVector.size(); ++j) {
        if (!isNA(datVector(j))) {
          if (datVector(j) > Type(0)) {
            true_count++;
          }
        }
      }

      // only compute observation likelihood if data vector has more than 1 non-zero element
      if (true_count > 1) {

        // generate index of non-zero values
        vector <int> datIdx(true_count);
        true_count = 0;
        for (int j = 0; j < datVector.size(); ++j) {
          if (!isNA(datVector(j))) {
            if (datVector(j) > Type(0)) {
              datIdx(true_count) = j;
              true_count++;
            }
          }
        }

        // extract non-zero elements
        vector <Type> datNonZero   = datVector(datIdx);
        vector <Type> alphaNonZero = alphaVector(datIdx);

        // negative log-likelihood
        nll += -ddirichlet(datNonZero,
                           alphaNonZero,
                           true);
      }
    }

    // report fitted variance-covariance matrix
    REPORT(MVNcov);

    // report random parameters
    REPORT(alphaMatrix);
    REPORT(probabilityMatrix);


  } else {
    error("Unknown model type");
  }
  return nll;
}
