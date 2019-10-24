void print_int_array(int *a, int num_elements){
	/*Print an array of integers*/
	int i;
	for (i = 0; i < num_elements; i++)
	{
		printf("%d\n", a[i]);
	}
	printf("\n");
}

void print_double_array(double *a, int num_elements){
	/*Print an array of doubles*/
	int i;
	for (i = 0; i < num_elements; i++)
	{
		printf("%.12f\n", a[i]);
	}
	printf("\n");
}

double distance_metric(int *data, int *simulation){
	/*Comupte the SMC distance metric between data and simulation

	Parameters
	----------------
	data : an array of length N_DATA, which is the data which we are using for
		Bayesian inference
	simulation : an array of length N_DATA, which is a simulated dataset

	Returns
	----------------
	distance : a double, the distance metric between the data and simulation
		arrays

	*/
	int i;
	int sum_data = 0;
	int sum_simulation = 0;
	double distance;

	for (i = 0; i < N_DATA; i++) {
		sum_data = sum_data + data[i];
		sum_simulation += simulation[i];
	}

	distance = (double)abs(sum_data - sum_simulation)/((double)N_DATA);
	return distance;
}

double kernel_pdf(double theta_old, double theta_new){
	/*The probability density of a perturbation kernel

	Parameters
	----------------
	theta_old : the value of the parameter at the previous time step
	theta_old : the value of the parameter at the current time step

	Returns
	----------------
	Transition probability from theta_old to theta_new

	*/
	return gsl_ran_gaussian_pdf(theta_new - theta_old, KERNEL_SD);
}

double prior_pdf(double theta){
	/*The probability density of a parameter under the prior*/
	return gsl_ran_beta_pdf(theta, PRIOR_ALPHA, PRIOR_BETA);
}

double update_distance_threshold(double *distance){
	/*Given a vector of distances, return a quantile defined by
	QUANTILE_ACCEPT_DISTANCE*/
	gsl_sort(distance, 1, N_PARTICLES);
	return gsl_stats_quantile_from_sorted_data(distance, 1, N_PARTICLES,
		QUANTILE_ACCEPT_DISTANCE);
}

int weighted_choice(gsl_rng *r, double *weight){
	/*Sample from an array of normalised weights

	Parameters
	----------------
	r : A GSL random number generator
	weight : An array of weights which sum to 1 (note: function will not check
		this for efficiency)

	Returns
	----------------
	An index from weights corresponding to a weighted sample.

	If function fails, returns -1
	*/

	double u = RND;
	int i;
	double up_to = 0.0;
	for (i = 0; i < N_PARTICLES; i++) {
		up_to += weight[i];
		if (up_to >= u) return i;
	}
	printf("Error in weighted_choice(). Weights not normalised!\n");
	return -1;
}

void write_particles_to_csv(double **theta_particle){
	/*Write a double array of dimensions (N_ROUNDS_SMC X N_PARTICLES) to the file
	OUTFILE_NAME*/

	FILE *outfile_pointer;
	outfile_pointer = fopen(OUTFILE_NAME, "w");

	int i, j;
	for (j = 0; j < N_PARTICLES; j++) {
		for (i = 0; i < N_ROUNDS_SMC; i++) {
			if (i < N_ROUNDS_SMC - 1) fprintf(outfile_pointer,"%.8f,", theta_particle[i][j]);
			else fprintf(outfile_pointer,"%.8f\n", theta_particle[i][j]);
		}
	}
	fclose(outfile_pointer);
}
