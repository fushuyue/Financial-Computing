#include <stdlib.h>
#include "MatrixStructures.h"
#include <random>
#include <iostream>
#include <time.h>

#define NUM_NON_ZERO_ENTRIES 10
#define NUM_ROWS 5
#define NUM_COLUMNS 5

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_spmatrix.h>

void test_that_student_has_boost_installed();
void test_that_student_has_GSL_installed();
void init_matrix(SparseMatrix & A);


int main(void)
{
	/*******Initial Checks*********/
	//test_that_student_has_boost_installed();
	//test_that_student_has_GSL_installed();

	SparseMatrix A, B, C;
	
	/*Initializing a sparse matrix*/
	srand(1); //uncomment durring development to get a consistant matrix.
	
	srand(time(NULL)); //I will use this to initialize matricies for grading.
	
	
	A._matrix = new SparseMatrix::element[NUM_NON_ZERO_ENTRIES];
	A._len_element_list = NUM_NON_ZERO_ENTRIES;
	A._m = NUM_ROWS;
	A._n = NUM_COLUMNS;
	
	B._matrix = new SparseMatrix::element[NUM_NON_ZERO_ENTRIES];
	B._len_element_list = NUM_NON_ZERO_ENTRIES;
	B._m = NUM_COLUMNS;
	B._n = NUM_ROWS;
	
	C._m = NUM_COLUMNS;
	C._n = NUM_ROWS;
	
	init_matrix(A);
	while(A.is_valid() == false)
	{
		//std::cout<<"init A is not valid, picking another\n";
		init_matrix(A);
	}
	
	init_matrix(B);
	while(B.is_valid() == false)
	{
		//std::cout<<"init B is not valid, picking another\n";
		init_matrix(B);
	}
	
	/* in the variable C, store the result of A*B */
	//A.print();
	C = A*B; //Do something to me!! Make me work, Please?
	/**********************************************/
	/**********************************************/

	std::cout<<"Print Matrix A: \n";
	A.print();
	std::cout<<"Print Matrix B: \n";
	B.print();
	std::cout<<"Checking if C is valid\n";
	if(C.is_valid())
		std::cout<<"C is valid\n";
	else
	{
		std::cout<<"C is NOT valid, exiting\n";
		return 1;
	}
	
	std::cout<<"Print Matrix C: \n";
	C.print();
	
};

void init_matrix(SparseMatrix & A)
{
	float tmp_value = 0.0;
	int tmp_i, tmp_j;
	for(int ii = 0; ii < A._len_element_list; ii++)
	{
		tmp_value = (rand()%1000 + 1)*0.1;
		tmp_i = rand()%A._m;
		tmp_j = rand()%A._n;
		
		A._matrix[ii]._value = tmp_value;
		A._matrix[ii]._i = tmp_i;
		A._matrix[ii]._j = tmp_j;
	}
};



void test_that_student_has_boost_installed()
{
	std::cout<<"From Boost test function\n";
	using namespace boost::numeric::ublas;
	mapped_matrix<double> m (3, 3, 3 * 3);
	for (unsigned i = 0; i < m.size1 (); ++ i)
		for (unsigned j = 0; j < m.size2 (); ++ j)
			m (i, j) = 3 * i + j;
	std::cout << m << std::endl;
};


void test_that_student_has_GSL_installed()
{
	std::cout<<"From GSL test function\n";
	
	gsl_spmatrix *A = gsl_spmatrix_alloc(5, 4); /* triplet format */
	gsl_spmatrix *C;
	size_t i, j;
	
	/* build the sparse matrix */
	gsl_spmatrix_set(A, 0, 2, 3.1);
	gsl_spmatrix_set(A, 0, 3, 4.6);
	gsl_spmatrix_set(A, 1, 0, 1.0);
	gsl_spmatrix_set(A, 1, 2, 7.2);
	gsl_spmatrix_set(A, 3, 0, 2.1);
	gsl_spmatrix_set(A, 3, 1, 2.9);
	gsl_spmatrix_set(A, 3, 3, 8.5);
	gsl_spmatrix_set(A, 4, 0, 4.1);
	
	printf("printing all matrix elements:\n");
	for (i = 0; i < 5; ++i)
		for (j = 0; j < 4; ++j)
			printf("A(%zu,%zu) = %g\n", i, j,
				   gsl_spmatrix_get(A, i, j));
	
	/* print out elements in triplet format */
	printf("matrix in triplet format (i,j,Aij):\n");
	for (i = 0; i < A->nz; ++i)
		printf("(%zu, %zu, %.1f)\n", A->i[i], A->p[i], A->data[i]);
	
	/* convert to compressed column format */
	C = gsl_spmatrix_compcol(A);
	
	printf("matrix in compressed column format:\n");
	printf("i = [ ");
	for (i = 0; i < C->nz; ++i)
		printf("%zu, ", C->i[i]);
	printf("]\n");
	
	printf("p = [ ");
	for (i = 0; i < C->size2 + 1; ++i)
		printf("%zu, ", C->p[i]);
	printf("]\n");
	
	printf("d = [ ");
	for (i = 0; i < C->nz; ++i)
		printf("%g, ", C->data[i]);
	printf("]\n");
	
	gsl_spmatrix_free(A);
	gsl_spmatrix_free(C);
}
