//based on https://github.com/stan-dev/example-models/blob/master/misc/multivariate-probit/probit-multi-good.stan
//all errors my own
data {
	// nY_all: number of individuals with all viruses tested
	int<lower=1> nY_all;
	// nY_SBPV: number of individuals with new plus SBPV	
	int<lower=1> nY_SBPV;
	// nY_reduced: number of individuals with new	
	int<lower=1> nY_reduced;
	// nV_all: number of viruses tested for
	int<lower=1> nV_all;
	// nV_SBPV: number of viruses tested for - no ABPV
	int<lower=1> nV_SBPV;
	// nV_reduced: number of viruses tested for - no ABPV, no SBPV
	int<lower=1> nV_reduced;
	// nL: number of locations
	int<lower=1> nL;
	// nS: number of species
	int<lower=1> nS;
	// nLS: number of species-location combinations
	int<lower=1> nLS;
	// nX: number of location level predictors
	int<lower=1> nX;
	// L_all: location membership for individuals with all viruses tested
	int<lower=1,upper=nL> L_all[nY_all];
	// L_SBPV: location membership for individuals with new plus SBPV
	int<lower=1,upper=nL> L_SBPV[nY_SBPV];
	// L_reduced: location membership for individuals with new	
	int<lower=1,upper=nL> L_reduced[nY_reduced];
	// S_all: location membership for individuals with all viruses tested
	int<lower=1,upper=nS> S_all[nY_all];
	// S_SBPV: location membership for individuals with new plus SBPV
	int<lower=1,upper=nS> S_SBPV[nY_SBPV];
	// S_reduced: location membership for individuals with new	
	int<lower=1,upper=nS> S_reduced[nY_reduced];
	// LS_all: location membership for individuals with all viruses tested
	int<lower=1,upper=nLS> LS_all[nY_all];
	// LS_SBPV: location membership for individuals with new plus SBPV
	int<lower=1,upper=nLS> LS_SBPV[nY_SBPV];
	// LS_reduced: location membership for individuals with new	
	int<lower=1,upper=nLS> LS_reduced[nY_reduced];
	// Y_all: outcome matrix
	int<lower=0,upper=1> Y_all[nY_all,nV_all];
	// Y_SBPV: outcome matrix
	int<lower=0,upper=1> Y_SBPV[nY_SBPV,nV_SBPV];
	// Y_reduced: outcome matrix
	int<lower=0,upper=1> Y_reduced[nY_reduced,nV_reduced];
	// X_all: predictor matrix
	vector[nX] X_all[nY_all];
	// X_SBPV: predictor matrix
	vector[nX] X_SBPV[nY_SBPV];
	// X_reduced: predictor matrix
	vector[nX] X_reduced[nY_reduced];
	
	// for prediction - counts
	int<lower=1> nY_pasc_all;
	int<lower=1> nY_pasc_SBPV;
	int<lower=1> nY_pasc_ABPV;
	int<lower=1> nY_luc_all;
	int<lower=1> nY_luc_SBPV;
	int<lower=1> nY_luc_ABPV;
	int<lower=1> nY_terr_all;
	int<lower=1> nY_terr_SBPV;
	int<lower=1> nY_terr_ABPV;
	int<lower=1> nY_other_all;
	int<lower=1> nY_other_SBPV;
	int<lower=1> nY_other_ABPV;
	// for prediction - covariates
	vector[nX] X_pred_pasc_all[nY_pasc_all];
	vector[nX] X_pred_pasc_SBPV[nY_pasc_SBPV];
	vector[nX] X_pred_pasc_ABPV[nY_pasc_ABPV];
	vector[nX] X_pred_luc_all[nY_luc_all];
	vector[nX] X_pred_luc_SBPV[nY_luc_SBPV];
	vector[nX] X_pred_luc_ABPV[nY_luc_ABPV];
	vector[nX] X_pred_terr_all[nY_terr_all];
	vector[nX] X_pred_terr_SBPV[nY_terr_SBPV];
	vector[nX] X_pred_terr_ABPV[nY_terr_ABPV];
	vector[nX] X_pred_other_all[nY_other_all];
	vector[nX] X_pred_other_SBPV[nY_other_SBPV];
	vector[nX] X_pred_other_ABPV[nY_other_ABPV];
	// for prediction - location
	int<lower=1,upper=nL> L_pred_pasc_all[nY_pasc_all];
	int<lower=1,upper=nL> L_pred_pasc_SBPV[nY_pasc_SBPV];
	int<lower=1,upper=nL> L_pred_pasc_APBV[nY_pasc_ABPV];
	int<lower=1,upper=nL> L_pred_luc_all[nY_luc_all];
	int<lower=1,upper=nL> L_pred_luc_SBPV[nY_luc_SBPV];
	int<lower=1,upper=nL> L_pred_luc_APBV[nY_luc_ABPV];
	int<lower=1,upper=nL> L_pred_terr_all[nY_terr_all];
	int<lower=1,upper=nL> L_pred_terr_SBPV[nY_terr_SBPV];
	int<lower=1,upper=nL> L_pred_terr_APBV[nY_terr_ABPV];
	int<lower=1,upper=nL> L_pred_other_all[nY_other_all];
	int<lower=1,upper=nL> L_pred_other_SBPV[nY_other_SBPV];
	int<lower=1,upper=nL> L_pred_other_APBV[nY_other_ABPV];
	// for prediction - species
	int<lower=1,upper=nS> S_pred_pasc_all[nY_pasc_all];
	int<lower=1,upper=nS> S_pred_pasc_SBPV[nY_pasc_SBPV];
	int<lower=1,upper=nS> S_pred_pasc_APBV[nY_pasc_ABPV];
	int<lower=1,upper=nS> S_pred_luc_all[nY_luc_all];
	int<lower=1,upper=nS> S_pred_luc_SBPV[nY_luc_SBPV];
	int<lower=1,upper=nS> S_pred_luc_APBV[nY_luc_ABPV];
	int<lower=1,upper=nS> S_pred_terr_all[nY_terr_all];
	int<lower=1,upper=nS> S_pred_terr_SBPV[nY_terr_SBPV];
	int<lower=1,upper=nS> S_pred_terr_APBV[nY_terr_ABPV];
	int<lower=1,upper=nS> S_pred_other_all[nY_other_all];
	int<lower=1,upper=nS> S_pred_other_SBPV[nY_other_SBPV];
	int<lower=1,upper=nS> S_pred_other_APBV[nY_other_ABPV];
	// for prediction - location/species
	int<lower=1,upper=nLS> LS_pred_pasc_all[nY_pasc_all];
	int<lower=1,upper=nLS> LS_pred_pasc_SBPV[nY_pasc_SBPV];
	int<lower=1,upper=nLS> LS_pred_pasc_APBV[nY_pasc_ABPV];
	int<lower=1,upper=nLS> LS_pred_luc_all[nY_luc_all];
	int<lower=1,upper=nLS> LS_pred_luc_SBPV[nY_luc_SBPV];
	int<lower=1,upper=nLS> LS_pred_luc_APBV[nY_luc_ABPV];
	int<lower=1,upper=nLS> LS_pred_terr_all[nY_terr_all];
	int<lower=1,upper=nLS> LS_pred_terr_SBPV[nY_terr_SBPV];
	int<lower=1,upper=nLS> LS_pred_terr_APBV[nY_terr_ABPV];
	int<lower=1> LS_pred_other_all[nY_other_all];
	int<lower=1,upper=nLS> LS_pred_other_SBPV[nY_other_SBPV];
	int<lower=1,upper=nLS> LS_pred_other_APBV[nY_other_ABPV];
}
parameters {
	// alpha: global latent variable intercept for each virus
	vector[nV_all] alpha;
	// sigmaL: scale parameters for locations
	vector<lower=0>[nV_all] sigmaL;
	// sigmaS: scale parameters for host species
	vector<lower=0>[nV_all] sigmaS;
	// sigmaLS: scale parameters for host species-location
	vector<lower=0>[nV_all] sigmaLS;
	// L_Omega: correlation of errors
	cholesky_factor_corr[nV_all] L_Omega;
	// beta: matrix of coefficients
	matrix[nV_all, nX] beta;
	
	//non-centered means
	// muL: specific parameters for location
	vector[nV_all] muL_tilde[nL];
	// muS: specific parameters for species
	vector[nV_all] muS_tilde[nS];
	// muLS: specific parameters for species-location
	vector[nV_all] muLS_tilde[nLS];
	
	//from git
	real<lower=0,upper=1> u_all[nY_all,nV_all]; // nuisance that absorbs inequality constraints
	real<lower=0,upper=1> u_SBPV[nY_SBPV,nV_SBPV]; // nuisance that absorbs inequality constraints
	real<lower=0,upper=1> u_reduced[nY_reduced,nV_reduced]; // nuisance that absorbs inequality constraints
}
transformed parameters {
	// muL: specific parameters for location
	vector[nV_all] muL[nL];
	// muS: specific parameters for species
	vector[nV_all] muS[nS];
	// muLS: specific parameters for species-location
	vector[nV_all] muLS[nLS];
	
	//calculated non-centered means
	for (l in 1:nL)
		for (n in 1:nV_all)
			muL[l,n] = muL_tilde[l,n] * sqrt(sigmaL[n]);
	for (s in 1:nS)
		for (n in 1:nV_all)
			muS[s,n] = muS_tilde[s,n] * sqrt(sigmaS[n]);
	for (m in 1:nLS)
		for (n in 1:nV_all)
			muLS[m,n] = muLS_tilde[m,n] * sqrt(sigmaLS[n]);
}
model {
	//priors
	//uniform prior over error covariate matrix
	L_Omega ~ lkj_corr_cholesky(1);
	//exponential prior over scale vectors
	sigmaL ~ exponential(0.5);
	sigmaS ~ exponential(0.5);
	sigmaLS ~ exponential(0.5);
	//diffuse normal prior for intercept
	alpha ~ logistic(0,1);
	//standard normal priors for beta coefficients - weakly regularising
	to_vector(beta) ~ normal(0,1);
	//shrinkage estimation of location values - assumed 0 correlation between viruses
	for (l in 1:nL)
		muL_tilde[l] ~ normal(0,1);
	//shrinkage estimation of species values - assumed 0 correlation between viruses
	for (s in 1:nS)
		muS_tilde[s] ~ normal(0,1);
	//shrinkage estimation of species-location values - assumed 0 correlation between viruses
	for (m in 1:nLS)
		muLS_tilde[m] ~ normal(0,1);
  	//implicit: u is iid standard uniform a priori
  	//model - splitting is possible due the rules of choleski decomposition 
  { // likelihood
    for (n in 1:nY_reduced) {
      vector[nV_reduced] mu;
      vector[nV_reduced] z;
      real prev;
      mu = head(beta * X_reduced[n], nV_reduced) + head(alpha, nV_reduced) + head(muS[S_reduced[n]], nV_reduced) + head(muL[L_reduced[n]], nV_reduced) + head(muLS[LS_reduced[n]], nV_reduced);
      prev = 0;
      for (d in 1:nV_reduced) { // Phi and inv_Phi may overflow and / or be numerically inaccurate
        real bound; // threshold at which utility = 0
        bound = Phi( -(mu[d] + prev) / L_Omega[d,d]  );
        if (Y_reduced[n,d] == 1) {
          real t;
          t = bound + (1 - bound) * u_reduced[n,d];
          z[d] = inv_Phi(t);       // implies utility is positive
          target += log1m(bound);  // Jacobian adjustment
        }
        else {
          real t;
          t = bound * u_reduced[n,d];
          z[d] = inv_Phi(t);     // implies utility is negative
          target += log(bound);  // Jacobian adjustment
        }
        if (d < nV_reduced) prev = L_Omega[d+1,1:d] * head(z, d);
        // Jacobian adjustments imply z is truncated standard normal
        // thus utility --- mu + L_Omega * z --- is truncated multivariate normal
      }
    }
  }
  { // likelihood
  for (n in 1:nY_all) {
      vector[nV_all] mu;
      vector[nV_all] z;
      real prev;
      mu = beta * X_all[n] + alpha + muS[S_all[n]] + muL[L_all[n]] + muLS[LS_all[n]];
      prev = 0;
      for (d in 1:nV_all) { // Phi and inv_Phi may overflow and / or be numerically inaccurate
        real bound; // threshold at which utility = 0
        bound = Phi( -(mu[d] + prev) / L_Omega[d,d]  );
        if (Y_all[n,d] == 1) {
          real t;
          t = bound + (1 - bound) * u_all[n,d];
          z[d] = inv_Phi(t);       // implies utility is positive
          target += log1m(bound);  // Jacobian adjustment
        }
        else {
          real t;
          t = bound * u_all[n,d];
          z[d] = inv_Phi(t);     // implies utility is negative
          target += log(bound);  // Jacobian adjustment
        }
        if (d < nV_all) prev = L_Omega[d+1,1:d] * head(z, d);
        // Jacobian adjustments imply z is truncated standard normal
        // thus utility --- mu + L_Omega * z --- is truncated multivariate normal
      }
    }
  }
  { // likelihood
    for (n in 1:nY_SBPV) {
      vector[nV_SBPV] mu;
      vector[nV_SBPV] z;
      real prev;
      mu = head(beta * X_SBPV[n], nV_SBPV) + head(alpha, nV_SBPV) + head(muS[S_SBPV[n]], nV_SBPV) + head(muL[L_SBPV[n]], nV_SBPV) + head(muLS[LS_SBPV[n]], nV_SBPV);
      prev = 0;
      for (d in 1:nV_SBPV) { // Phi and inv_Phi may overflow and / or be numerically inaccurate
        real bound; // threshold at which utility = 0
        bound = Phi( -(mu[d] + prev) / L_Omega[d,d]  );
        if (Y_SBPV[n,d] == 1) {
          real t;
          t = bound + (1 - bound) * u_SBPV[n,d];
          z[d] = inv_Phi(t);       // implies utility is positive
          target += log1m(bound);  // Jacobian adjustment
        }
        else {
          real t;
          t = bound * u_SBPV[n,d];
          z[d] = inv_Phi(t);     // implies utility is negative
          target += log(bound);  // Jacobian adjustment
        }
        if (d < nV_SBPV) prev = L_Omega[d+1,1:d] * head(z, d);
        // Jacobian adjustments imply z is truncated standard normal
        // thus utility --- mu + L_Omega * z --- is truncated multivariate normal
      }
    }
  }
}
generated quantities {
  //correlation matrix
  corr_matrix[nV_all] Omega;
  //simulation variables
  int<lower=0> y_pred_pasc_all[6];
  int<lower=0> y_pred_terr_all[6];
  int<lower=0> y_pred_luc_all[6];
  int<lower=0> y_pred_other_all[6];
  int<lower=0> y_pred_pasc_ABPV;
  int<lower=0> y_pred_terr_ABPV;
  int<lower=0> y_pred_luc_ABPV;
  int<lower=0> y_pred_other_ABPV;
  int<lower=0> y_pred_pasc_SBPV[2];
  int<lower=0> y_pred_terr_SBPV[2];
  int<lower=0> y_pred_luc_SBPV[2];
  int<lower=0> y_pred_other_SBPV[2];
  //final count variables
  int<lower=0> pasc_RLV;
  int<lower=0> pasc_LMV;
  int<lower=0> pasc_MV1;
  int<lower=0> pasc_MV2;
  int<lower=0> pasc_SBPV;
  int<lower=0> pasc_ABPV;
  int<lower=0> terr_RLV;
  int<lower=0> terr_LMV;
  int<lower=0> terr_MV1;
  int<lower=0> terr_MV2;
  int<lower=0> terr_SBPV;
  int<lower=0> terr_ABPV;
  int<lower=0> luc_RLV;
  int<lower=0> luc_LMV;
  int<lower=0> luc_MV1;
  int<lower=0> luc_MV2;
  int<lower=0> luc_SBPV;
  int<lower=0> luc_ABPV;
  int<lower=0> other_RLV;
  int<lower=0> other_LMV;
  int<lower=0> other_MV1;
  int<lower=0> other_MV2;
  int<lower=0> other_SBPV;
  int<lower=0> other_ABPV;
  
  //define zeros for multi_normal_rng
  vector[6] zeros;
  zeros = rep_vector(0,6);
  
  //calculate correlation matrix
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  
  //initialise counting variables
  for (i in 1:6) {
    y_pred_pasc_all[i] = 0;
    y_pred_terr_all[i] = 0;
    y_pred_luc_all[i] = 0;
    y_pred_other_all[i] = 0;
  }
  for (i in 1:2) {
    y_pred_pasc_SBPV[i] = 0;
    y_pred_terr_SBPV[i] = 0;
    y_pred_luc_SBPV[i] = 0;
    y_pred_other_SBPV[i] = 0;
  }
  y_pred_pasc_ABPV = 0;
  y_pred_terr_ABPV = 0;
  y_pred_luc_ABPV = 0;
  y_pred_other_ABPV = 0;
  
  
  //generate predictions for individuals in pools but not tested in main dataset
  //and calculate sums - we ignore the known values, statistically less satisfying
  //but algebraically simpler
  
  //pascuorum
  for (n in 1:nY_pasc_all) {
  	vector[6] mu;
  	vector[6] sims;
  	mu = beta * X_pred_pasc_all[n] + alpha + muL[L_pred_pasc_all[n]] + muS[S_pred_pasc_all[n]] + muLS[LS_pred_pasc_all[n]];
  	sims = multi_normal_rng(mu, Omega);
  	for (i in 1:6) {
  	  if (sims[i] >= 0) {
  	    y_pred_pasc_all[i] = y_pred_pasc_all[i] + 1;
  	  }
  	}
  }
  for (n in 1:nY_pasc_SBPV) {
  	vector[6] mu;
  	vector[6] sims;
  	mu = beta * X_pred_pasc_SBPV[n] + alpha + muL[L_pred_pasc_SBPV[n]] + muS[S_pred_pasc_SBPV[n]] + muLS[LS_pred_pasc_SBPV[n]];
  	sims = multi_normal_rng(mu, Omega);
  	for (i in 5:6) {
  	  if (sims[i] >= 0) {
  	    y_pred_pasc_SBPV[(i-4)] = y_pred_pasc_SBPV[(i-4)] + 1;
  	  }
  	}
  }
  for (n in 1:nY_pasc_ABPV) {
  	vector[6] mu;
  	vector[6] sims;
  	mu = beta * X_pred_pasc_ABPV[n] + alpha + muL[L_pred_pasc_APBV[n]] + muS[S_pred_pasc_APBV[n]] + muLS[LS_pred_pasc_APBV[n]];
  	sims = multi_normal_rng(mu, Omega);
  	if (sims[6] >= 0) {
  	  y_pred_pasc_ABPV = y_pred_pasc_ABPV + 1;
  	}
  }
  //terrestris
  for (n in 1:nY_terr_all) {
  	vector[6] mu;
  	vector[6] sims;
  	mu = beta * X_pred_terr_all[n] + alpha + muL[L_pred_terr_all[n]] + muS[S_pred_terr_all[n]] + muLS[LS_pred_terr_all[n]];
  	sims = multi_normal_rng(mu, Omega);
  	for (i in 1:6) {
  	  if (sims[i] >= 0) {
  	    y_pred_terr_all[i] = y_pred_terr_all[i] + 1;
  	  }
  	}
  }
  for (n in 1:nY_terr_SBPV) {
  	vector[6] mu;
  	vector[6] sims;
  	mu = beta * X_pred_terr_SBPV[n] + alpha + muL[L_pred_terr_SBPV[n]] + muS[S_pred_terr_SBPV[n]] + muLS[LS_pred_terr_SBPV[n]];
  	sims = multi_normal_rng(mu, Omega);
  	for (i in 5:6) {
  	  if (sims[i] >= 0) {
  	    y_pred_terr_SBPV[(i-4)] = y_pred_terr_SBPV[(i-4)] + 1;
  	  }
  	}
  }
  for (n in 1:nY_terr_ABPV) {
  	vector[6] mu;
  	vector[6] sims;
  	mu = beta * X_pred_terr_ABPV[n] + alpha + muL[L_pred_terr_APBV[n]] + muS[S_pred_terr_APBV[n]] + muLS[LS_pred_terr_APBV[n]];
  	sims = multi_normal_rng(mu, Omega);
  	if (sims[6] >= 0) {
  	  y_pred_terr_ABPV = y_pred_terr_ABPV + 1;
  	}
  }
  //lucorum
  for (n in 1:nY_luc_all) {
  	vector[6] mu;
  	vector[6] sims;
  	mu = beta * X_pred_luc_all[n] + alpha + muL[L_pred_luc_all[n]] + muS[S_pred_luc_all[n]] + muLS[LS_pred_luc_all[n]];
  	sims = multi_normal_rng(mu, Omega);
  	for (i in 1:6) {
  	  if (sims[i] >= 0) {
  	    y_pred_luc_all[i] = y_pred_luc_all[i] + 1;
  	  }
  	}
  }
  for (n in 1:nY_luc_SBPV) {
  	vector[6] mu;
  	vector[6] sims;
  	mu = beta * X_pred_luc_SBPV[n] + alpha + muL[L_pred_luc_SBPV[n]] + muS[S_pred_luc_SBPV[n]] + muLS[LS_pred_luc_SBPV[n]];
  	sims = multi_normal_rng(mu, Omega);
  	for (i in 5:6) {
  	  if (sims[i] >= 0) {
  	    y_pred_luc_SBPV[(i-4)] = y_pred_luc_SBPV[(i-4)] + 1;
  	  }
  	}
  }
  for (n in 1:nY_luc_ABPV) {
  	vector[6] mu;
  	vector[6] sims;
  	mu = beta * X_pred_luc_ABPV[n] + alpha + muL[L_pred_luc_APBV[n]] + muS[S_pred_luc_APBV[n]] + muLS[LS_pred_luc_APBV[n]];
  	sims = multi_normal_rng(mu, Omega);
  	if (sims[6] >= 0) {
  	  y_pred_luc_ABPV = y_pred_luc_ABPV + 1;
  	}
  }
  //others
  for (n in 1:nY_other_all) {
  	vector[6] mu;
  	vector[6] sims;
  	if (LS_pred_other_all[n] > nLS) {
  	  mu = beta * X_pred_other_all[n] + alpha + muL[L_pred_other_all[n]] + muS[S_pred_other_all[n]] + multi_normal_rng(zeros, diag_matrix(sigmaLS));
  	} else {
  	  mu = beta * X_pred_other_all[n] + alpha + muL[L_pred_other_all[n]] + muS[S_pred_other_all[n]] + muLS[LS_pred_other_all[n]];
  	}
  	sims = multi_normal_rng(mu, Omega);
  	for (i in 1:6) {
  	  if (sims[i] >= 0) {
  	    y_pred_other_all[i] = y_pred_other_all[i] + 1;
  	  }
  	}
  }
  for (n in 1:nY_other_SBPV) {
  	vector[6] mu;
  	vector[6] sims;
  	mu = beta * X_pred_other_SBPV[n] + alpha + muL[L_pred_other_SBPV[n]] + muS[S_pred_other_SBPV[n]] + muLS[LS_pred_other_SBPV[n]];
  	sims = multi_normal_rng(mu, Omega);
  	for (i in 5:6) {
  	  if (sims[i] >= 0) {
  	    y_pred_other_SBPV[(i-4)] = y_pred_other_SBPV[(i-4)] + 1;
  	  }
  	}
  }
  for (n in 1:nY_other_ABPV) {
  	vector[6] mu;
  	vector[6] sims;
  	mu = beta * X_pred_other_ABPV[n] + alpha + muL[L_pred_other_APBV[n]] + muS[S_pred_other_APBV[n]] + muLS[LS_pred_other_APBV[n]];
  	sims = multi_normal_rng(mu, Omega);
  	if (sims[6] >= 0) {
  	  y_pred_other_ABPV = y_pred_other_ABPV + 1;
  	}
  }
  
  //sum values for output
  pasc_RLV = y_pred_pasc_all[1];
  pasc_LMV = y_pred_pasc_all[2];
  pasc_MV1 = y_pred_pasc_all[3];
  pasc_MV2 = y_pred_pasc_all[4];
  pasc_SBPV = y_pred_pasc_all[5] + y_pred_pasc_SBPV[1];
  pasc_ABPV = y_pred_pasc_all[6] + y_pred_pasc_SBPV[2] + y_pred_pasc_ABPV;
  terr_RLV = y_pred_terr_all[1];
  terr_LMV = y_pred_terr_all[2];
  terr_MV1 = y_pred_terr_all[3];
  terr_MV2 = y_pred_terr_all[4];
  terr_SBPV = y_pred_terr_all[5] + y_pred_terr_SBPV[1];
  terr_ABPV = y_pred_terr_all[6] + y_pred_terr_SBPV[2] + y_pred_terr_ABPV;
  luc_RLV = y_pred_luc_all[1];
  luc_LMV = y_pred_luc_all[2];
  luc_MV1 = y_pred_luc_all[3];
  luc_MV2 = y_pred_luc_all[4];
  luc_SBPV = y_pred_luc_all[5] + y_pred_luc_SBPV[1];
  luc_ABPV = y_pred_luc_all[6] + y_pred_luc_SBPV[2] + y_pred_luc_ABPV;
  other_RLV = y_pred_other_all[1];
  other_LMV = y_pred_other_all[2];
  other_MV1 = y_pred_other_all[3];
  other_MV2 = y_pred_other_all[4];
  other_SBPV = y_pred_other_all[5] + y_pred_other_SBPV[1];
  other_ABPV = y_pred_other_all[6] + y_pred_other_SBPV[2] + y_pred_other_ABPV;
}
