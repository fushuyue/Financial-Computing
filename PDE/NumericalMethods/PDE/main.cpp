#include <iostream>
#include <fstream>

#include "do_nm.h"

int main(void)
{
	std::ofstream output_file;
	output_file.open("out.csv");
	int N = 20;
	int M = 10;
	double r = 0.02;
	double sigma = .2;
	double K = 20;
	double S_max = 100;
	double T = 20;
	//double ** do_euler_explicit(int N, int M, float r, double sigma, double K, double S_max, double T)
	double ** result_table = NULL;
	result_table = do_euler_explicit(N, M, r, sigma, K, S_max, T);
	std::cout<<"The Result!\n";
	for(int i =0; i < M; i++)
	{
		for(int j = 0; j<N-1; j++)
		{
	  std::cout<<result_table[i][j]<<" ";
	  output_file<<result_table[i][j]<<',';
		}
		std::cout<<result_table[i][N-1]<<std::endl;
		output_file<<result_table[i][N-1]<<"\n";
	}
	
}