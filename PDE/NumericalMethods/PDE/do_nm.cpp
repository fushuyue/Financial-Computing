#include <iostream>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <cmath>
#include <gsl/gsl_linalg.h>

double max(double a,double b)
{
	if(a>b)
		return a;
	else
		return b;
}
double min(double a,double b)
{
	if(a>b)
		return b;
	else
		return a;
}



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
	double xbar = log(S_max);
	double delta_t = T/M;
	double h = xbar*2.0/(N);
	double beta = r-0.5*sigma*sigma;
	double k1 = delta_t*sigma*sigma/(2*h*h)-beta*delta_t/(2.0*h);
	double k2 = 1-delta_t*sigma*sigma/(h*h)-delta_t*r;
	double k3 = delta_t*sigma*sigma/(2*h*h)+beta*delta_t/(2.0*h);
	
	gsl_matrix * t_m = gsl_matrix_alloc(N,N);
	gsl_matrix_set_zero(t_m);
	for(int i = 0; i < N; i++)
	{
		
		gsl_matrix_set(t_m, i, i, k2);
	}
	for(int i = 0; i < N-1; i++)
	{
		
		gsl_matrix_set(t_m, i+1, i, k1);
	}
	for(int i = 0; i < N-1; i++)
		gsl_matrix_set(t_m, i, i+1, k3);
	
	// print
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
	for(int i=0;i<N;i++)
	{
		double temp;
		double sPriceTemp = -xbar+i*h;
		if((exp(sPriceTemp)-K)>0)
			temp=exp(sPriceTemp)-K;
		else
			temp=0;
		gsl_vector_set(s_v_old, i, temp);
		
	}
	
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
		{
			/* boundary condion two: if price is too large, use the put-call parity */
			double optionPrice = gsl_vector_get(s_v_new, i);
			
			if (optionPrice<0.00001)
				gsl_vector_set(s_v_new, i, 0);
			
			if (optionPrice>exp(-xbar+i*h))
				gsl_vector_set(s_v_new, i, exp(-xbar+i*h));
			
			if (i==(N-1))
				gsl_vector_set(s_v_new, i,(exp(-xbar+i*h)-K*exp(-r*(t_j*delta_t))) );
			
			
			std::cout<<gsl_vector_get(s_v_new, i)<<" ";
			
		}
		std::cout<<std::endl;
		
		memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
		
		gsl_vector_swap(s_v_new, s_v_old);
	}
	
	std::cout<<std::endl;
	std::cout<<std::endl;
	
	std::cout<<"/**************************************************************************************************************/"<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	gsl_matrix_free(t_m);
	gsl_vector_free(s_v_old);
	gsl_vector_free(s_v_new);
	
	/***** Above is an example, replace with whatever you want to use to build this method ************************/
	/**************************************************************************************************************/
	
	return result_table;
};

double ** do_euler_implicitet(int N, int M, float r, double sigma, double K, double S_max, double T)
{
	double ** result_table;
	result_table = new double *[M];
	for(int i =0; i<M; i++)
	{
		result_table[i] = new double[N];
		for(int j = 0; j<N; j++)
			result_table[i][j] = 0;
	}
	
	// Build Example Matrix
	double xbar = log(S_max);
	double delta_t = T/M;
	double h = xbar*2/(N);
	double beta = r-0.5*sigma*sigma;
	double k1 = -delta_t*sigma*sigma/(2*h*h)+beta*delta_t/(2*h);
	double k2 = 1+delta_t*sigma*sigma/(h*h)+delta_t*r;
	double k3 = -delta_t*sigma*sigma/(2*h*h)-beta*delta_t/(2*h);
	
	gsl_matrix * t_m = gsl_matrix_alloc(N,N);
	gsl_matrix_set_zero(t_m);
	for(int i = 0; i < N; i++)
	{
		
		gsl_matrix_set(t_m, i, i, k2);
	}
	for(int i = 0; i < N-1; i++)
	{
		
		gsl_matrix_set(t_m, i+1, i, k1);
	}
	for(int i = 0; i < N-1; i++)
		gsl_matrix_set(t_m, i, i+1, k3);
	
	// print
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
	for(int i=0;i<N;i++)
	{
		double temp;
		double sPriceTemp = -xbar+i*h;
		if((exp(sPriceTemp)-K)>0)
			temp=exp(sPriceTemp)-K;
		else
			temp=0;
		gsl_vector_set(s_v_old, i, temp);
		
	}
	
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
	
	int s;
	gsl_permutation * p = gsl_permutation_alloc (N);
	gsl_linalg_LU_decomp (t_m, p, &s);
	
	// this is the loop that steps through time
	for(int t_j = 0; t_j < M; t_j++)
	{
		
		gsl_linalg_LU_solve (t_m, p, s_v_old, s_v_new);
		//gsl_permutation_free (p);
		
		for(int i = 0; i < N; i++)
		{
			/* boundary condion two: if price is too large, use the put-call parity */
			double optionPrice = gsl_vector_get(s_v_new, i);
			if (optionPrice<0.00001)
				gsl_vector_set(s_v_new, i, 0);
			if (optionPrice>exp(-xbar+i*h))
				gsl_vector_set(s_v_new, i, exp(-xbar+i*h));
			if (i==(N-1))
				gsl_vector_set(s_v_new, i,(exp(-xbar+i*h)-K*exp(-r*(t_j*delta_t))) );
			
			std::cout<<gsl_vector_get(s_v_new, i)<<" ";
			
		}
		std::cout<<std::endl;
		
		memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
		
		gsl_vector_swap(s_v_new, s_v_old);
	}
	
	std::cout<<std::endl;
	std::cout<<std::endl;
	
	std::cout<<"/**************************************************************************************************************/"<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	gsl_matrix_free(t_m);
	gsl_vector_free(s_v_old);
	gsl_vector_free(s_v_new);
	
	/***** Above is an example, replace with whatever you want to use to build this method ************************/
	/**************************************************************************************************************/
	
	return result_table;
	
};

double ** do_crank_nicolson(int N, int M, float r, double sigma, double K, double S_max, double T)
{
	double ** result_table;
	result_table = new double *[M];
	for(int i =0; i<M; i++)
	{
		result_table[i] = new double[N];
		for(int j = 0; j<N; j++)
			result_table[i][j] = 0;
	}
	
	// Build Example Matrix
	double xbar = log(S_max);
	double delta_t = T/M;
	double h = xbar*2/(N);
	double beta = r-0.5*sigma*sigma;
	double k1 = -sigma*sigma/(4*h*h)+beta/(4*h);
	double k2 = 1.0/delta_t+sigma*sigma/(2*h*h)+r/2;
	double k3 = -sigma*sigma/(4*h*h)-beta/(4*h);
	double k1star = sigma*sigma/(4*h*h)-beta/(4*h);
	double k2star = 1.0/delta_t-sigma*sigma/(2*h*h)-r/2;
	double k3star = sigma*sigma/(4*h*h)+beta/(4*h);
	gsl_matrix * t_m = gsl_matrix_alloc(N,N);
	gsl_matrix_set_zero(t_m);
	gsl_matrix * t_n = gsl_matrix_alloc(N,N);
	gsl_matrix_set_zero(t_n);
	
	for(int i = 0; i < N; i++)
	{
		
		gsl_matrix_set(t_m, i, i, k2);
		gsl_matrix_set(t_n, i, i, k2star);
	}
	for(int i = 0; i < N-1; i++)
	{
		gsl_matrix_set(t_n, i+1, i, k1star);
		gsl_matrix_set(t_m, i+1, i, k1);
	}
	for(int i = 0; i < N-1; i++)
	{
		gsl_matrix_set(t_n, i, i+1, k3star);
		gsl_matrix_set(t_m, i, i+1, k3);
	}
	// print
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
			std::cout<<gsl_matrix_get(t_m, i, j)<<" ";
		std::cout<<std::endl;
	}
	
	// Solution Vector
	gsl_vector * s_v_old = gsl_vector_alloc(N);
	gsl_vector * s_v_new = gsl_vector_alloc(N);
	gsl_vector * s_v_temp = gsl_vector_alloc(N);
	
	gsl_vector_set_zero(s_v_old);
	gsl_vector_set_zero(s_v_new);
	gsl_vector_set_zero(s_v_temp);
	
	//example initial value
	
	for(int i=0;i<N;i++)
	{
		double temp;
		double sPriceTemp = -xbar+i*h;
		if((exp(sPriceTemp)-K)>0)
			temp=exp(sPriceTemp)-K;
		else
			temp=0;
		gsl_vector_set(s_v_old, i, temp);
	}
	
	//gsl_vector_set(s_v_old, 0, k2star*d[0]+k3star*d[1]);
	//gsl_vector_set(s_v_old, N-1,k1star*d[N-2]+k2star*d[N-1]+k3star*(exp(xbar)-K*exp(-r*T)));
	
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
	
	int s=0;
	gsl_matrix*t_m_inverse=gsl_matrix_alloc(N,N);
	gsl_matrix_set_zero(t_m_inverse);
	gsl_permutation * p = gsl_permutation_alloc (N);
	gsl_linalg_LU_decomp (t_m, p, &s);
	gsl_linalg_LU_invert (t_m, p, t_m_inverse);
	
	// this is the loop that steps through time
	for(int t_j = 0; t_j < M; t_j++)
	{
		gsl_blas_dgemv(CblasNoTrans, 1.0, t_n, s_v_old, 0.0, s_v_temp);
		gsl_vector_set(s_v_temp, N-1,gsl_vector_get(s_v_temp, N-1)+k3star*(exp(xbar)-K*exp(-r*T)));
		
		gsl_blas_dgemv(CblasNoTrans, 1.0, t_m_inverse, s_v_temp, 0.0, s_v_new);
		//gsl_linalg_LU_solve (t_m, p, s_v_temp, s_v_new);
		
		
		//gsl_permutation_free (p);
		
		std::cout<<"new:"<<std::endl;
		for(int i = 0; i < N; i++)
		{
			/* boundary condion two: if price is too large, use the put-call parity */
			double optionPrice = gsl_vector_get(s_v_new, i);
			std::cout<<gsl_vector_get(s_v_new, i)<<" ";
			
			if (optionPrice<0.00001)
				gsl_vector_set(s_v_new, i, 0);
			//if (optionPrice>exp(-xbar+i*h))
			//	gsl_vector_set(s_v_new, i, exp(-xbar+i*h));
			if (i==(N-1))
				gsl_vector_set(s_v_new, i,(exp(-xbar+i*h)-K*exp(-r*(t_j*delta_t))) );
			
			std::cout<<gsl_vector_get(s_v_new, i)<<" ";
			
		}
		
		
		memcpy(result_table[t_j], s_v_new->data, s_v_new->size*sizeof(double));
		gsl_vector_swap(s_v_new, s_v_old);
	}
	
	std::cout<<std::endl;
	std::cout<<std::endl;
	
	std::cout<<"/**************************************************************************************************************/"<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	std::cout<<std::endl;
	gsl_matrix_free(t_m);
	gsl_vector_free(s_v_old);
	gsl_vector_free(s_v_new);
	
	/***** Above is an example, replace with whatever you want to use to build this method ************************/
	/**************************************************************************************************************/
	
	return result_table;
};
