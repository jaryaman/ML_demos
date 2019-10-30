/*
Performing approximate Bayesian computation sequential Monte Carlo (Toni et al.
2009) for linear regression.

Synthetic data is generated by `ground_truth_and_analysis.ipynb`.

This script writes particles.csv, where each row corresponds to a particle and
each column corresponds to a round of SMC.

Defining #DEBUG_MODE will silence all writing to stdout. One may then add
printf statements in the code, and perhaps write the output to file as:
`./run.sh > output.txt`

Parameter ordering convention:
0 - gradient
1 - intercept
2 - standard deviation

Author: Juvid Aryaman
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fit.h>


#define N_DATA 30
#define N_PARAMETERS 3

#define N_PARTICLES 2000
#define N_ROUNDS_SMC 25
#define QUANTILE_ACCEPT_DISTANCE 0.8

#define RND gsl_rng_uniform(r)
#define SEED 1
#define DISTANCE_THRESHOLD_INIT_GRADIENT 50
#define DISTANCE_THRESHOLD_INIT_INTERCEPT 50
#define DISTANCE_THRESHOLD_INIT_SIGMA 50

#define X_DATA_FILENAME "x.csv"
#define Y_DATA_FILENAME "y.csv"

//#define DEBUG_MODE

#include "smc.h"
#include "lin_reg.h"

int main(int argc, char *argv[]) {

/* set up GSL RNG */
gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
/* end of GSL setup */

gsl_rng_set(r, SEED);

/////////////////////////
/*Read data*/
/////////////////////////

FILE *data_pointer_x, *data_pointer_y;

data_pointer_x = fopen(X_DATA_FILENAME, "r");
data_pointer_y = fopen(Y_DATA_FILENAME, "r");

double data_x[N_DATA];
double data_y[N_DATA];
int i, j, read_error_status_x, read_error_status_y;
for (i=0; i < N_DATA; i++){
	read_error_status_x = fscanf(data_pointer_x, "%lf\n", &data_x[i]);
	read_error_status_y = fscanf(data_pointer_y, "%lf\n", &data_y[i]);
}
if (read_error_status_x != 1){printf("Error reading X data\n"); return 0;}
if (read_error_status_y != 1){printf("Error reading Y data\n"); return 0;}

/////////////////////////
/*Initialise variables*/
/////////////////////////

/*Fit a linear model to the data, which will be used as summary statistics of
the data*/
double gradient_fit_data, intercept_fit_data, sigma_fit_data, cov00, cov01, cov11, sumsq;
int gsl_fit_return_value;
gsl_fit_return_value = gsl_fit_linear(data_x, 1, data_y, 1, N_DATA,
																			&intercept_fit_data, &gradient_fit_data,
																			&cov00, &cov01, &cov11, &sumsq);
if (gsl_fit_return_value != 0) {printf("Fit failed.\n"); return -1;}
sigma_fit_data = sqrt(sumsq/(N_DATA-2));

#ifndef DEBUG_MODE
	printf("gradient ML = %.8f\n", gradient_fit_data);
	printf("intercept ML = %.8f\n", intercept_fit_data);
	printf("sigma ML = %.8f\n", sigma_fit_data);
#endif



/*Make a (N_PARAMETERS X N_ROUNDS_SMC X N_PARTICLES) array to store all
particles at all rounds of SMC*/
double*** theta_particle;
theta_particle = (double***) malloc(N_PARAMETERS * sizeof(double**));
for (i = 0; i < N_PARAMETERS; i++){
	theta_particle[i] = (double**) malloc(N_ROUNDS_SMC * sizeof(double*));
}
for (i = 0; i < N_PARAMETERS; i++){
	for (j = 0; j < N_ROUNDS_SMC; j++){
		theta_particle[i][j] = (double*) malloc(N_PARTICLES * sizeof(double));
	}
}

double** distance;
distance = (double**) malloc(N_PARAMETERS * sizeof(double*));
for (i = 0; i < N_PARAMETERS; i++) {
	distance[i] = (double*) malloc(N_PARTICLES * sizeof(double));
}
double **distance_threshold_all = malloc(N_PARAMETERS * sizeof(double*));
for (i = 0; i < N_PARAMETERS; i++) {
	distance_threshold_all[i] = (double*) malloc(N_ROUNDS_SMC * sizeof(double));
}

double distance_threshold[] = {DISTANCE_THRESHOLD_INIT_GRADIENT,
															 DISTANCE_THRESHOLD_INIT_INTERCEPT,
														   DISTANCE_THRESHOLD_INIT_SIGMA};

double *simulated_data = (double*) malloc(N_DATA * sizeof(double));
double *weight = malloc(N_PARTICLES * sizeof(double));

double weight_normalizer = 0.0;

int time_smc=0; // an index of each round of SMC
int param_index_chosen, prior_violated;
int particle_index;
/////////////////////////
/*Perform ABC SMC*/
/////////////////////////

/*For every round of SMC*/
for (time_smc = 0; time_smc < N_ROUNDS_SMC; time_smc++) {
	#ifndef DEBUG_MODE
		printf("Round %d of SMC\n", time_smc);
	#endif
	for (i = 0; i < N_PARAMETERS; i++) {
		distance_threshold_all[i][time_smc] = distance_threshold[i];
	}

	/*Draw or perturb a particle and compute distance*/
	for (particle_index = 0; particle_index < N_PARTICLES; particle_index++) {
		// printf("%d\n", particle_index);
		for (i = 0; i < N_PARAMETERS; i++) {
			// reset distance of particle along each dimension
			distance[i][particle_index] = distance_threshold[i] + 1.0;
		}
		while ((distance[0][particle_index] > distance_threshold[0])||
					 (distance[1][particle_index] > distance_threshold[1])||
				   (distance[2][particle_index] > distance_threshold[2])) {
			if (time_smc == 0) {
				// Sample from the prior
				sample_prior(r, theta_particle, particle_index);
			}
			else{
				/*Sample from the old weights*/
				param_index_chosen = weighted_choice(r, weight);
				if ((param_index_chosen < 0)||(param_index_chosen >= N_PARTICLES)) {
					printf("Error in param_index_chosen\n");
					printf("time_smc = %d\n", time_smc);
					return -1;
				}

				perturb_particle(r, theta_particle, time_smc, param_index_chosen,
												 particle_index);

				// Check if prior support is 0
				prior_violated = check_prior_violated(theta_particle, time_smc,
																							particle_index);
				if(prior_violated == 1) continue;
				}

			simulate_dataset(r, theta_particle, data_x, simulated_data, time_smc,
											 particle_index);




			// Compute distance between data and simulation
			distance_metric_sum_stats(simulated_data,
															 data_x,
															 gradient_fit_data,
															 intercept_fit_data,
															 sigma_fit_data,
														   distance, particle_index);
			//distance[particle_index] = distance_metric_sum_res(simulated_data, data_y);
			// distance[particle_index] = distance_metric_sum_abs_res(simulated_data, data_y);

		}
	}

	#ifndef DEBUG_MODE
		printf("Particles sampled.\n");
	#endif


	/*Compute weights*/
	if (time_smc==0){ for (i = 0; i < N_PARTICLES; i++) weight[i] = 1.0;}
	else{
		weight_normalizer = 0.0;
		for (particle_index = 0; particle_index < N_PARTICLES; particle_index++) {
			weight_normalizer +=  weight[particle_index]*kernel_pdf();
		}
		// print_double_array(weight, N_PARTICLES);
		// printf("\n" );
		// printf("weight_normalizer=%f\n", weight_normalizer);
		// printf("kernel_pdf()=%f\n", kernel_pdf());
		for (particle_index = 0; particle_index < N_PARTICLES; particle_index++) {
			weight[particle_index] = prior_pdf(theta_particle, time_smc,
																				 particle_index)/weight_normalizer;
		}

	}

	/*Normalise weights*/
	weight_normalizer = 0.0;
	for (i = 0; i < N_PARTICLES; i++)	weight_normalizer += weight[i];
	for (i = 0; i < N_PARTICLES; i++)	weight[i] = weight[i]/weight_normalizer;

	/* Resample weights*/
	distance_threshold[0] = update_distance_threshold(distance[0]);
	distance_threshold[1] = update_distance_threshold(distance[1]);
	distance_threshold[2] = update_distance_threshold(distance[2]);


}

#ifndef DEBUG_MODE
	printf("Writing particles to file\n");
#endif
	write_particles_to_csv(theta_particle);

	char *dist_filename = "distances.txt";
	write_2d_double_array_to_csv(distance_threshold_all, N_PARAMETERS, N_ROUNDS_SMC, dist_filename);
#ifndef DEBUG_MODE
	printf("Done!\n");
#endif

return 0; //return from main
} //close main
