data {
     int<lower=0> N;			              // number of observations
     int<lower=1> P;                    // number of fixed effects
     int<lower=0> I;                    // number of patients
     int<lower=1> nr; 		   	          // number patients random effects
     int<lower=0> Time;		   	          // number of visits
     int<lower=0> M;                    // number of treatments
     int<lower=1,upper=M> treat[N];	    // treatment id
     int<lower=1,upper=Time> time[N];  	// visit id
     int<lower=1,upper=I> patient[N];   // patient id
     row_vector[P] X[N];                // fixed effects design matrix
     row_vector[nr] Z[N];               // patient random effect design matrix
     vector[N] visual;      		        // visual acuity
     vector[N] visual0;	   	            // visual acuity in time 0
}

parameters {
	   vector[P] beta;    		// fixed coefficients
	   vector[nr] zr[I];		      // patient intercepts
	   real<lower=0> sigma_e; // residual std
	   real<lower=0> sigma_p;	// patient std
}

model {
      // priors
      for (k in 1:I)
        zr[k] ~ normal(0,sigma_p);
      // likelihood
      for (i in 1:N)
        visual[i] ~ normal(X[i] * beta + 
        Z[i] * zr[patient[i]],sigma_e);
}      
      	  