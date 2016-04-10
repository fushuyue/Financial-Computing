#include <iostream>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

double ** do_euler_explicit(int N, int M, float r, double sigma, double K, double S_max, double T)
{
	
	/**************************************************************************************************************/
	/***** Below is an example, replace with whatever you want to use to build this method ************************/
	
	double ** result_table;
	result_table = new double *[M];
	for(int i =0; i<M; i++)
	{
		result_table[i] = new double[N];
		for(int j = 0; j<N; j++)
			result_table[i][j] = 0;
	}
	
	// Build Example Matrix
	gsl_matrix * t_m = gsl_matrix_alloc(N,N);
	gsl_matrix_set_zero(t_m);
	for(int i = 0; i < N; i++)
		gsl_matrix_set(t_m, i, i, 0.25);
	
	for(int i = 0; i < N-1; i++)
		gsl_matrix_set(t_m, i+1, i, -0.25);
	
	for(int i = 0; i < N-1; i++)
		gsl_matrix_set(t_m, i, i+1, -0.25);
	
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
			std::cout<<gsl_matrix_get(t_m, i, j)<<" ";
		std::cout<<std::endl;
	}
	
	// Solution Vector
	gsl_vector * s_v_old = gsl_vector_alloc(N);
	gsl_vector * s_v_new = gsl_vector_alloc(N);
	
	gsl_vector_set_zero(s_v_old);
	gsl_vector_set_zero(s_v_new);
	
	//example initial value
	gsl_vector_set(s_v_old, N/2, 2.5);
	gsl_vector_set(s_v_old, N/2 - 1, 1.5);
	gsl_vector_set(s_v_old, N/2 + 1, 1.5);
	
	std::cout<<"New\n";
	for(int i = 0; i < N; i++)
		std::cout<<gsl_vector_get(s_v_new, i)<<" ";
	std::cout<<std::endl;
	
	std::cout<<"Old\n";
	for(int i = 0; i < N; i++)
		std::cout<<gsl_vector_get(s_v_old, i)<<" ";
	std::cout<<std::endl;
	
	// printing out the initial condition
	std::cout<<"Starting\n";
	for(int i = 0; i < N; i++)
		std::cout<<gsl_vector_get(s_v_old, i)<<" ";
	std::cout<<std::endl;
	
	// this is the loop that steps through time
	for(int t_j = 0; t_j < M; t_j++)
	{
		// This is what actually does the matrix vector multiplication.
		gsl_blas_dgemv(CblasNoTrans, 1.0, t_m, s_v_old, 0.0, s_v_new);
		
		for(int i = 0; i < N; i++)
			std::cout<<gsl_vector_get(s_v_new, i)<<" ";
		std::cout<<std::endl;
		
		memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
		
		gsl_vector_swap(s_v_new, s_v_old);
	}
	
	gsl_matrix_free(t_m);
	gsl_vector_free(s_v_old);
	gsl_vector_free(s_v_new);
	
	/***** Above is an example, replace with whatever you want to use to build this method ************************/
	/**************************************************************************************************************/
	
	return result_table;
};

double ** do_euler_implicitet(int N, int M, float r, double sigma, double K, double S_max, double T)
{
	/**************************************************************************************************************/
	/***** Below is an example, replace with whatever you want to use to build this method ************************/
	
	double ** result_table;
	result_table = new double *[M];
	for(int i =0; i<M; i++)
	{
		result_table[i] = new double[N];
		for(int j = 0; j<N; j++)
			result_table[i][j] = 0;
	}
	
	/***** Above is an example, replace with whatever you want to use to build this method ************************/
	/**************************************************************************************************************/
	
	return result_table;
};

double ** do_crank_nicolson(int N, int M, float r, double sigma, double K, double S_max, double T)
{
	/**************************************************************************************************************/
	/***** Below is an example, replace with whatever you want to use to build this method ************************/
	
	double ** result_table;
	result_table = new double *[M];
	for(int i =0; i<M; i++)
	{
		result_table[i] = new double[N];
		for(int j = 0; j<N; j++)
			result_table[i][j] = 0;
	}
	
	/***** Above is an example, replace with whatever you want to use to build this method ************************/
	/**************************************************************************************************************/
	return result_table;
};