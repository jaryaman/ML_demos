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


	// double update_distance_threshold(double *distance){
	// 	/*Given a vector of distances, return a quantile defined by
	// 	QUANTILE_ACCEPT_DISTANCE*/
	// 	gsl_sort(distance, 1, N_PARTICLES);
	// 	return gsl_stats_quantile_from_sorted_data(distance, 1, N_PARTICLES,
	// 		QUANTILE_ACCEPT_DISTANCE);
	// }

void write_particles_to_csv(double ***theta_particle){
	/*Write a double array of dimensions (N_ROUNDS_SMC X N_PARTICLES) to the file
	OUTFILE_NAME*/

	FILE *outfile_pointer;

	int i, j, k;
	char *outfile_name = (char*)malloc(50 * sizeof(char));

	for (k = 0; k < N_PARAMETERS; k++) {
		sprintf(outfile_name, "particle_%d.csv", k);
		outfile_pointer = fopen(outfile_name, "w");
		for (j = 0; j < N_PARTICLES; j++) {
			for (i = 0; i < N_ROUNDS_SMC; i++) {
				if (i < N_ROUNDS_SMC - 1) fprintf(outfile_pointer,"%.8f,", theta_particle[k][i][j]);
				else fprintf(outfile_pointer,"%.8f\n", theta_particle[k][i][j]);
			}
		}
		fclose(outfile_pointer);
	}
}

void write_double_array_to_csv(double *arr, char *filename){
	/*Write a double array of length N_ELEMENTS to file*/

	FILE *outfile_pointer;
	int i;
	int num_elements = sizeof(arr) / sizeof(double);
	outfile_pointer = fopen(filename, "w");
	for (i = 0; i < num_elements; i++) {
		fprintf(outfile_pointer,"%.8f\n", arr[i]);
	}
	fclose(outfile_pointer);
}


void write_2d_double_array_to_csv(double **arr, int N_ROWS, int N_COLS, char *filename){
	/*Write a 2D double array of length N_ELEMENTS to file*/

	FILE *outfile_pointer;
	int i, j;

	outfile_pointer = fopen(filename, "w");
	for (i = 0; i < N_ROWS; i++) {
		for (j = 0; j < N_COLS; j++) {
			if (j < N_COLS-1) fprintf(outfile_pointer,"%.8f,", arr[i][j]);
			else fprintf(outfile_pointer,"%.8f\n", arr[i][j]);
		}
	}
	fclose(outfile_pointer);
}
