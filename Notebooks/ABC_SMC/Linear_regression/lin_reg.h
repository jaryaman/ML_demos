#define PRIOR_GRADIENT_LOWER 0.0
#define PRIOR_INTERCEPT_LOWER 3.0
#define PRIOR_SIGMA_LOWER 0.0

#define PRIOR_GRADIENT_UPPER 10.0
#define PRIOR_INTERCEPT_UPPER 500.0
#define PRIOR_SIGMA_UPPER 10.0

#define KERNEL_SD_GRADIENT 0.05
#define KERNEL_SD_INTERCEPT 5.0
#define KERNEL_SD_SIGMA 0.1

double unif_neg_pos(gsl_rng *r){
  /*Return a unif(-1,1)*/
  return 2.0*gsl_rng_uniform(r) - 1.0;
}


void sample_prior(gsl_rng *r, double ***theta_particle, int particle_index){
  /*Sample from prior for linear regression

  Parameters
	----------------
	r : A GSL random number generator
	theta_particle : An array of dimensions (N_PARAMETERS X N_ROUNDS_SMC X
    N_PARTICLES) containing parameter values at SMC time points (i.e. particles)
  particle_index : Index of a particle

	Returns
	----------------
	Augments theta_particle to sample from the prior for parameter, for particle
  particle_index, and adds this to the 0th round of SMC slot
  */
  theta_particle[0][0][particle_index] = (PRIOR_GRADIENT_UPPER -
    PRIOR_GRADIENT_LOWER)*gsl_rng_uniform(r) + PRIOR_GRADIENT_LOWER;
  theta_particle[1][0][particle_index] = (PRIOR_INTERCEPT_UPPER -
    PRIOR_INTERCEPT_LOWER)*gsl_rng_uniform(r) + PRIOR_INTERCEPT_LOWER;
  theta_particle[2][0][particle_index] = (PRIOR_SIGMA_UPPER -
    PRIOR_SIGMA_LOWER)*gsl_rng_uniform(r) + PRIOR_SIGMA_LOWER;
}

int check_prior_violated(double ***theta_particle, int time_smc,
                         int particle_index){
  /*Check if the support of the prior for any parameter, for a particular
  particle at time point time_smc, is 0

  Parameters
  ----------------
  theta_particle : an array of dimensions (N_PARAMETERS X N_ROUNDS_SMC X
    N_PARTICLES) containing parameter values at SMC time points (i.e. particles)
  time_smc : The time point in SMC
  particle_index : Index of a particle

  Returns
  ----------------
  1 if priors are violated, 0 otherwise

  */
  if ((theta_particle[0][time_smc][particle_index] < PRIOR_GRADIENT_LOWER) ||
      (theta_particle[0][time_smc][particle_index] > PRIOR_GRADIENT_UPPER)){
        return 1;
      }
  if ((theta_particle[1][time_smc][particle_index] < PRIOR_INTERCEPT_LOWER) ||
      (theta_particle[1][time_smc][particle_index] > PRIOR_INTERCEPT_UPPER)){
        return 1;
      }
  if ((theta_particle[2][time_smc][particle_index] < PRIOR_SIGMA_LOWER) ||
      (theta_particle[2][time_smc][particle_index] > PRIOR_SIGMA_UPPER)) {
        return 1;
      }
  return 0;
}


void perturb_particle(gsl_rng *r, double ***theta_particle, int time_smc,
                       int param_index_chosen, int particle_index){
  /* Perturb particles at a given SMC time point

  Parameters
  ----------------
  r : A GSL random number generator
  theta_particle : An array of dimensions (N_PARAMETERS X N_ROUNDS_SMC X
    N_PARTICLES) containing parameter values at SMC time points (i.e. particles)
  time_smc : The current time point for which a perturbed parametrization is to
    be produced
  param_index_chosen : The index of the parameter from (time_smc-1) to be
    perturbed
  particle_index : Index of a particle

  Returns
  ----------------
  Augments theta_particle at time point time_smc to contain perturbed parameters
  from time_smc-1

  */
  int i;
  double perturbation_kernel[] = {KERNEL_SD_GRADIENT,
                                  KERNEL_SD_INTERCEPT,
                                  KERNEL_SD_SIGMA};
  double u;
  for (i = 0; i < N_PARAMETERS; i++) {
    u = unif_neg_pos(r); // Unif(-1,1)
    theta_particle[i][time_smc][particle_index] =
      theta_particle[i][time_smc-1][param_index_chosen] +
      perturbation_kernel[i]*u;
  }
}


void simulate_dataset(gsl_rng *r, double ***theta_particle, double *data_x,
  double *simulated_data, int time_smc, int particle_index){
  /*Simulate a linear regression dataset and add to simulated_data

  Parameters
  ----------------
  r : A GSL random number generator
  theta_particle : an array of dimensions (N_PARAMETERS X N_ROUNDS_SMC X
    N_PARTICLES) containing parameter values at SMC time points (i.e. particles)
  data_x : an array corresponding to the independent variable x
  time_smc : The time point in SMC
  simulated_data : an array of length N_DATA, where each element is a regression
  against x, using parameters from theta_particle
  time_smc : The time point in SMC
  particle_index : Index of a particle

  Returns
  ----------------
  Augments simulated_data, filling it with a simulated dataset corresponding to
  particle `particle_index`, at `time_smc`
  */

  int i;
  double gradient, intercept, sigma;

  gradient = theta_particle[0][time_smc][particle_index];
  intercept = theta_particle[1][time_smc][particle_index];
  sigma = theta_particle[2][time_smc][particle_index];

  for (i = 0; i < N_DATA; i++) {
    simulated_data[i] = gradient*data_x[i] + intercept +
                        gsl_ran_gaussian(r, sigma);
  }
}


double distance_metric_sum_stats(double *simulated_data, double *data_x,
                       double gradient_fit_data, double intercept_fit_data,
                       double sigma_fit_data){
  /* Compute a distance metric between the data and the simulation as the sum
  of relative absolute distances between maximum-likelihood estimates of the
  three parameters of linear regression.

  Parameters
  ----------------
  simulated_data : An array of length N_DATA of simulated data
  data_x : An array of length N_DATA of the independent variable
  gradient_fit_data : ML fit of the gradient to the data
  intercept_fit_data : ML fit of the intercept to the data
  sigma_fit_data : ML fit of the standard deviation to the data


  Returns
  ----------------
  distance metric between the data and the simulation

  */

  double gradient_fit_sim, intercept_fit_sim, sigma_fit_sim, cov00, cov01, cov11, sumsq;
  int gsl_fit_return_value;
  double distance_metric;

  gsl_fit_return_value = gsl_fit_linear(data_x, 1, simulated_data, 1, N_DATA,
                                        &intercept_fit_sim, &gradient_fit_sim,
                                         &cov00, &cov01, &cov11, &sumsq);
  if (gsl_fit_return_value != 0) {printf("Fit failed.\n"); exit(99);}
  sigma_fit_sim = sqrt(sumsq/(N_DATA-2));
  distance_metric = fabs(gradient_fit_sim - gradient_fit_data)/gradient_fit_data +
                   fabs(intercept_fit_sim - intercept_fit_data)/intercept_fit_data +
                   fabs(sigma_fit_sim - sigma_fit_data)/sigma_fit_data;
  if (distance_metric < 0) {printf("Negative distance!\n");  exit(99);}

  return distance_metric;
}

double distance_metric_sum_sq_res(double *simulated_data, double *data_y){
  /* Compute a distance metric between the data and the simulation as the sum
  of squared residuals/N_DATA.

  NOTE: This is not a good distance metric for SMC because it will attempt to
  find a maximum-likelihood estimate for the gradient and intercept, which will
  cause the noise parameter to overfit.

  Parameters
  ----------------
  simulated_data : An array of length N_DATA of simulated data
  data_y : An array of length N_DATA of the dependent variable


  Returns
  ----------------
  distance metric between the data and the simulation

  */

  int i;
  double res = 0.0;
  for (i = 0; i < N_DATA; i++) {
    res += (data_y[i] - simulated_data[i])*(data_y[i] - simulated_data[i]);
  }
  return res/N_DATA;
}

double distance_metric_sum_abs_res(double *simulated_data, double *data_y){
  /* Compute a distance metric between the data and the simulation as the sum
  of absolute residuals.

  NOTE: This is not a good distance metric for SMC because it will attempt to
  find a maximum-likelihood estimate for the gradient and intercept, which will
  cause the noise parameter to overfit.

  Parameters
  ----------------
  simulated_data : An array of length N_DATA of simulated data
  data_y : An array of length N_DATA of the dependent variable


  Returns
  ----------------
  distance metric between the data and the simulation

  */

  int i;
  double res = 0.0;
  for (i = 0; i < N_DATA; i++) {
    res += fabs(data_y[i] - simulated_data[i]);
  }
  return res/N_DATA;
}

double kernel_pdf(){
	/*The probability density of a new parameter given an old parameter under the
	perturbation kernel

	Parameters
	----------------
	theta_old : the value of the parameter at the previous time step
	theta_old : the value of the parameter at the current time step

	Returns
	----------------
	Transition probability density from theta_old to theta_new

	*/

	return 1.0/(2.0*KERNEL_SD_GRADIENT)/(2.0*KERNEL_SD_INTERCEPT)/(2.0*KERNEL_SD_SIGMA);
}

double prior_pdf(double ***theta_particle, int time_smc, int particle_index){
	/*The probability density of a parameter under the prior

  Parameters
  ----------------
  theta_particle : an array of dimensions (N_PARAMETERS X N_ROUNDS_SMC X
    N_PARTICLES) containing parameter values at SMC time points (i.e. particles)
  time_smc : The time point in SMC
  particle_index : Index of a particle

  Returns
	----------------
  Prior probability of the parameter at time_smc for particle particle_index
  */
  int prior_is_violated;
  double prior = 1.0;
  prior_is_violated = check_prior_violated(theta_particle, time_smc,
                                           particle_index);
  if(prior_is_violated == 1){
    printf("Prior violated\n");
    return 0.0;
  }
  else{
    prior = prior/(PRIOR_GRADIENT_UPPER-PRIOR_GRADIENT_LOWER);
    prior = prior/(PRIOR_INTERCEPT_UPPER-PRIOR_INTERCEPT_LOWER);
    prior = prior/(PRIOR_SIGMA_UPPER-PRIOR_SIGMA_LOWER);
    return prior;
  }
}
