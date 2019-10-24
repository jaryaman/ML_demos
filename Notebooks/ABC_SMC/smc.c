#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics.h>

#define N_DATA 50
#define N_TRUTH 10

#define PRIOR_ALPHA 0.5
#define PRIOR_BETA 0.5

#define N_PARTICLES 5000
#define N_ROUNDS_SMC 50
#define KERNEL_SD 0.05
#define QUANTILE_ACCEPT_DISTANCE 0.8

#define RND gsl_rng_uniform(r)
#define SEED 1
#define DISTANCE_THRESHOLD_INIT 10

#define OUTFILE_NAME "particles.csv"

//#define DEBUG_MODE

#include "smc.h"

int main(int argc, char *argv[]) {


/* set up GSL RNG */
gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
/* end of GSL setup */

gsl_rng_set(r, SEED);

///////////////////////////////////////////////////////////////////
// READ DATA
///////////////////////////////////////////////////////////////////

FILE *data_pointer;

data_pointer = fopen("binom_data.csv", "r");

int data[N_DATA];
int i, j, read_error_status;

/*Make a (N_ROUNDS_SMC X N_PARTICLES) array to store all particles at all rounds
of SMC*/
double** theta_particle;
theta_particle = (double**) malloc(N_ROUNDS_SMC * sizeof(double*));
for (i = 0; i < N_ROUNDS_SMC; i++){
	theta_particle[i] = (double*) malloc(N_PARTICLES * sizeof(double));
}


for (i=0; i < N_DATA; i++){
	read_error_status = fscanf(data_pointer, "%d\n", &data[i]);
}
if (read_error_status != 1){printf("Error reading data\n"); return 0;}

///////////////////////////////////////////////////////////////////
// Perform ABC SMC
///////////////////////////////////////////////////////////////////

int time_smc=0; // an index of each round of SMC

double distance_threshold = DISTANCE_THRESHOLD_INIT;
int *simulated_data = (int*) malloc(N_DATA * sizeof(int));
double *distance = malloc(N_PARTICLES * sizeof(double));
double *weight = malloc(N_PARTICLES * sizeof(double));
double weight_normalizer = 0.0;

int param_index_chosen;

/*For every round of SMC*/
for (time_smc = 0; time_smc < N_ROUNDS_SMC; time_smc++) {
	#ifndef DEBUG_MODE
		printf("Round %d of SMC\n", time_smc);
	#endif

	/*Draw or perturb a particle and compute distance*/
	for (i = 0; i < N_PARTICLES; i++) {
		distance[i] = distance_threshold + 1.0; // reset distance of particle i
		while (distance[i] > distance_threshold) {
			if (time_smc == 0) {
				// Sample from the prior
				theta_particle[0][i] = gsl_ran_beta(r, PRIOR_ALPHA, PRIOR_BETA);
			}
			else{
				/*Sample from the old weights and perturb*/
				param_index_chosen = weighted_choice(r, weight);
				if ((param_index_chosen < 0)||(param_index_chosen >= N_PARTICLES)) {
					printf("Error in param_index_chosen\n"); return -1;
			}

				theta_particle[time_smc][i] =
					theta_particle[time_smc-1][param_index_chosen] +
					gsl_ran_gaussian(r, KERNEL_SD);

				// check if bounds of prior exceeded
				if((theta_particle[time_smc][i]<0) || (theta_particle[time_smc][i] > 1)) continue;
			}


			// Simulate a candidate dataset
			for (j = 0; j < N_DATA; j++) {
				simulated_data[j] = gsl_ran_binomial(r, theta_particle[time_smc][i], N_TRUTH);
			}

			distance[i] = distance_metric(data, simulated_data);

		}
	}
	#ifndef DEBUG_MODE
		printf("Particles sampled.\n");
	#endif

	/*Compute weights*/
	if (time_smc==0){ for (i = 0; i < N_PARTICLES; i++) weight[i] = 1.0;}
	else{
		weight_normalizer = 0.0;
		for (i = 0; i < N_PARTICLES; i++) {
			weight_normalizer += weight[i]*kernel_pdf(theta_particle[time_smc-1][i],
				theta_particle[time_smc][i]);
		}

		for (i = 0; i < N_PARTICLES; i++) {
			weight[i] = prior_pdf(theta_particle[time_smc][i])/weight_normalizer;
		}

	}


	/*Normalise weights*/
	weight_normalizer = 0.0;
	for (i = 0; i < N_PARTICLES; i++)	weight_normalizer += weight[i];
	for (i = 0; i < N_PARTICLES; i++)	weight[i] = weight[i]/weight_normalizer;


	/* Resample weights*/
	distance_threshold = update_distance_threshold(distance);

}

// // printf("\n");
//print_double_array(weight, N_PARTICLES);

#ifndef DEBUG_MODE
	printf("Writing particles to file\n");
#endif
	write_particles_to_csv(theta_particle);
#ifndef DEBUG_MODE
	printf("Done!\n");
#endif

return 0; //return from main
} //close main
