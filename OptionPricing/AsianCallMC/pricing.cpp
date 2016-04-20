//
//  pricing.cpp
//  IE525Project
//
//  Created by fushuyue on 4/19/16.
//  Copyright © 2016 Monte Carlo Simulation. All rights reserved.
//
#include "interface.h"
#include <iostream>
#include <cmath>

double callOptionPrice(double S,double t,double X,double r,double sigma,double q)
{
	
	double d1=(log(S/X) + (r-q+sigma*sigma/2)*t)/(sigma*sqrt(t));
	double d2=(log(S/X) + (r-q-sigma*sigma/2)*t)/(sigma*sqrt(t));
	return Normalcdf(d1)*S*exp(-q*t)- Normalcdf(d2)*X*exp(-r*t);
}

double control_monte_carlo_call_price(const int& n, const double& S, const double& K, const double& r, const double& q, const double& v, const double& T, double m)
{
	double time_interval=(double)T/m;
	double sigma_hat_square=(2*m+1)/(3*m)*v*v;
	double q_hat=q+0.5*v*v-0.5*sigma_hat_square;
	double T_hat=0.5*(T+time_interval);
	double x_bar=0;double y_bar=0;
	double memory1[1000],memory2[1000];
	
	for (int i=0;i<1000;i++)
	{
		double sum1=0,sum2=1;
		double S_adjust = S;//*exp(time_interval*(r-0.5*v*v-q));
		//cout<<S_adjust<<endl;
		for (int j=0;j<m;j++)
		{
			//cout<<random<<endl;
			S_adjust=S_adjust*exp(sqrt(v*v*time_interval)*gaussian_box_muller())*exp(time_interval*(r-0.5*v*v-q));
			sum1=S_adjust+sum1;
			sum2=sum2*S_adjust;
		}
		// memory1 is the arithemtic average price y
		// memory2 is the geometric average price x
		
		memory1[i]=sum1/m;
		memory2[i]=pow(sum2,1.0/m);
		
		memory1[i]=max(memory1[i]-K,0);
		memory2[i]=max(memory2[i]-K,0);
		
		x_bar=x_bar+memory2[i];
		y_bar=y_bar+memory1[i];
	}
	
	x_bar = x_bar/1000;
	y_bar=y_bar/1000;
	//std::cout << "The call price without control variate is :"<<y_bar<<std::endl;
	double top=0;
	double bottom=0;
	
	for (int i=0;i<1000;i++)
	{
		//std::cout<<memory2[i]<<" ";
		top=top+(memory2[i]-x_bar)*(memory1[i]-y_bar);
		bottom=bottom+pow(memory2[i]-x_bar,2);
	}
	//std::cout << "botton :"<<bottom<<std::endl;
	double b=top/bottom;
	std::cout<<"b= "<<b<<std::endl;
	
	double E_X=exp(r*T_hat)*callOptionPrice(S,T_hat,K,r,sqrt(sigma_hat_square),q_hat);
	//std::cout<<"E_X= "<<E_X<<"\n";
	double se_x=0,se_y=0;
	double x_final=0;
	double y_final=0;
	for (int i=0;i<n;i++)
	{
		double S_adjust=S;
		double sum1=0,sum2=1;
		for (int j=0;j<m;j++)
		{
			double random=gaussian_box_muller();
			S_adjust=S_adjust*exp(v*sqrt(time_interval)*random)*exp(time_interval*(r-0.5*v*v-q));
			sum1=S_adjust+sum1;
			sum2=sum2*S_adjust;
		}
		double result1=max((double)(sum1/m-K),0.0);
		double result2=max(pow(sum2,1.0/m)-K,0.0);
		
		double Y=result1+b*(E_X-result2);
		double X=result1;
		y_final=Y+y_final;
		x_final=X+x_final;
		se_y=Y*Y+se_y;
		se_x=X*X+se_x;
	}
	double se1,se2;
	se1=sqrt((se_y-n*(y_final/n)*(y_final/n))/(n-1))/sqrt(n);
	se2=sqrt((se_x-n*(x_final/n)*(x_final/n))/(n-1))/sqrt(n);
	
	std::cout << "The call price without control variate is :"<<exp(-r*T)*x_final/n<<std::endl;
	std::cout << "The call price with control variate is :"<<exp(-r*T)*y_final/n<<std::endl;
	std::cout << "The Standard Error without control variate is :"<<se2<<std::endl;
	std::cout << "The Standard Error control variate is :"<<se1<<std::endl;
	std::cout << "The confidence interval of 95% is : ["<<exp(-r*T)*y_final-1.96*se2<<" , "<<exp(-r*T)*y_final+1.96*se2<<"]"<<std::endl;
	std::cout << "The range of confidence interval is : "<<2*1.96*se2<<std::endl;
	
	return exp(-r*T)*y_final;
}


double Asian_call_qmc(const int& n, const double& S, const double& K, const double& r, const double& q, const double& v, const double& T, double m,double l){
	
	int D = m;
	char p[20] = "file.txt";
	double time_interval=(double)T/m;
	double sigma_hat_square=(2*m+1)/(3*m)*v*v;
	double q_hat=q+0.5*v*v-0.5*sigma_hat_square;
	double T_hat=0.5*(T+time_interval);
	double x_bar=0;double y_bar=0;
	double se_x_l=0,se_y_l=0;
	double x_final_l=0;
	double y_final_l=0;
	for(int ll=0;ll<l;ll++)
	{
		double **P = sobol_points(n,D,p);
		double *memory1,*memory2;
		memory1=new double [1000];
		memory2=new double [1000];
	
		for (int i=0;i<1000;i++)
		{
			double sum1=0,sum2=1;
			double S_adjust = S;
			for (int j=0;j<m;j++)
			{
				double random = Uniform()+P[i][j];
				random = fmod(random,1);
				random = inverse(random);
				S_adjust=S_adjust*exp(sqrt(v*v*time_interval)*random)*exp(time_interval*(r-0.5*v*v-q));
				sum1=S_adjust+sum1;
				sum2=sum2*S_adjust;
			}
			// memory1 is the arithemtic average price y
			// memory2 is the geometric average price x
			
			memory1[i]=sum1/m;
			memory2[i]=pow(sum2,1.0/m);
			
			memory1[i]=max(memory1[i]-K,0);
			memory2[i]=max(memory2[i]-K,0);
			
			x_bar=x_bar+memory2[i];
			y_bar=y_bar+memory1[i];
		}
		
		x_bar = x_bar/1000;
		y_bar=y_bar/1000;
		
		double top=0;
		double bottom=0;
		
		for (int i=0;i<1000;i++)
		{
			//std::cout<<memory2[i]<<" ";
			top=top+(memory2[i]-x_bar)*(memory1[i]-y_bar);
			bottom=bottom+pow(memory2[i]-x_bar,2);
		}
		//std::cout << "botton :"<<bottom<<std::endl;
		double b=top/bottom;
		std::cout<<"b= "<<b<<std::endl;
		
		double E_X=exp(r*T_hat)*callOptionPrice(S,T_hat,K,r,sqrt(sigma_hat_square),q_hat);
		//std::cout<<"E_X= "<<E_X<<"\n";
		double se_x=0,se_y=0;
		double x_final=0;
		double y_final=0;
		for (int i=0;i<n;i++)
		{
			double S_adjust=S;
			double sum1=0,sum2=1;
			for (int j=0;j<m;j++)
			{
				double random = Uniform()+P[i][j];
				random = fmod(random,1);
				random = inverse(random);
				S_adjust=S_adjust*exp(sqrt(v*v*time_interval)*random)*exp(time_interval*(r-0.5*v*v-q));
				sum1=S_adjust+sum1;
				sum2=sum2*S_adjust;
			}
			double result1=max((double)(sum1/m-K),0.0);
			double result2=max(pow(sum2,1.0/m)-K,0.0);
			
			double Y=result1+b*(E_X-result2);
			double X=result1;
			y_final=Y+y_final;
			x_final=X+x_final;
			se_y=Y*Y+se_y;
			se_x=X*X+se_x;
		}
		double se1,se2;
		se1=sqrt((se_y-n*(y_final/n)*(y_final/n))/(n-1))/sqrt(n);
		se2=sqrt((se_x-n*(x_final/n)*(x_final/n))/(n-1))/sqrt(n);
		x_final=exp(-r*T)*x_final/n;y_final=exp(-r*T)*y_final/n;
		x_final_l+=x_final;
		y_final_l+=y_final;
	}
	
	std::cout << "The call price without control variate is :"<<x_final_l/l<<std::endl;
	std::cout << "The call price with control variate is :"<<y_final_l/l<<std::endl;
	//std::cout << "The Standard Error without control variate is :"<<se2<<std::endl;
	//std::cout << "The Standard Error control variate is :"<<se1<<std::endl;
	//std::cout << "The confidence interval of 95% is : ["<<exp(-r*T)*y_final-1.96*se2<<" , "<<exp(-r*T)*y_final+1.96*se2<<"]"<<std::endl;
	//std::cout << "The range of confidence interval is : "<<2*1.96*se2<<std::endl;

	return 0;
}
