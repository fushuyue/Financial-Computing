//
//  main.cpp
//  hw3
//
//  Created by fushuyue on 2/25/16.
//  Copyright © 2016 Monte Carlo Simulation. All rights reserved.
//

#include <iostream>
#include <cmath>
const double k=100;
const double t=(double)1/52;
const double s0=100;
const double q=0;
const double sigma=0.3;
const double r=0.005;
const double discont=exp(-r*t);


double get_uniform(){
	return (double)(random()/(((double)pow(2,31)-1)));
}

double ar(){
	double u1,u2,u3,x,temp;
	do
	{
		u1=get_uniform();
		u2=get_uniform();
		x=-log(u1);
		temp= -0.5*(x-1)*(x-1);
	}while (u2>exp(temp));
	
	u3=get_uniform();
	if(u3<=0.5) return -x;
	else return x;

}

double get_st()
{
	double temp;
	temp=s0*exp(t*(r-q-0.5*0.3*0.3)+0.3*sqrt(t)*ar());
	
	return (double)temp;
}

double max(double a, double b){
	if(a>=b)
		return a;
	else
		return b;
}

double get_price()
{
	double st=get_st();
	if(st<=k) return discont*(-st+k);
	else return 0;
}

int main(int argc, const char * argv[]){
	
	double start_s=clock();
	long n=100000;
	// number of times for simulation
	//sscanf(argv[1],"%ld",&n);
	double *p=new double[n];
	for(long i=0;i<n;i++)
		p[i]=get_price();
	double sum=0;
	double sum1=0;
	for(long i=0;i<n;i++)
	{
		sum=sum+p[i];
		sum1=sum1+p[i]*p[i];
	}
	
	double se;
	se=sqrt((sum1-n*(sum/n)*(sum/n))/(n-1))/sqrt(n);
	std::cout<<n<<"	"<<sum/n<<"	"<<sum/n-101.65<<"	"<<se;
	double stop_s=clock();
	double time=(stop_s-start_s)/double(CLOCKS_PER_SEC);
	std::cout<<"	"<<time<<std::endl;
	
}
