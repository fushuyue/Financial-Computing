//
//  main.cpp
//  Shuyue_Fu_Assignment_10
//
//  Created by fushuyue on 11/22/15.
//  Copyright © 2015 Financial Computing. All rights reserved.
//
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

double up_factor, uptick_prob, downtick_prob, notick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;
// global arrays serve as memoization tools
double momi1[10000][10000];
double momi2[10000][10000];

// max function written by professor
double max(double a, double b) {
	return (b < a )? a:b;
}


double american_call_option(int k, int i, double current_stock_price)
{	double temp;
	if (momi1[k][i+no_of_divisions]!=-1) return momi1[k][i+no_of_divisions];
	if (k == no_of_divisions)
		return max(0.0, (current_stock_price - strike_price));
	else
	{
		// if the array is empty, put the computed value into the array
		if(momi1[k+1][i+no_of_divisions]==-1)
		{
			momi1[k+1][i+no_of_divisions]=american_call_option(k+1, i,current_stock_price);
		}
		if(momi1[k+1][i+1+no_of_divisions]==-1)
		{
			momi1[k+1][i+1+no_of_divisions]=american_call_option(k+1, i+1,current_stock_price*up_factor);
		}
		if(momi1[k+1][i-1+no_of_divisions]==-1)
		{
			momi1[k+1][i-1+no_of_divisions]=american_call_option(k+1, i-1,current_stock_price/up_factor);
		}
		
		// normal way of computing an American call
		temp=(uptick_prob*momi1[k+1][i+1+no_of_divisions]+downtick_prob*momi1[k+1][i-1+no_of_divisions]+ notick_prob*momi1[k+1][i+no_of_divisions])/R;
		momi1[k][i+no_of_divisions]=max(temp, (current_stock_price - strike_price));
		return momi1[k][i+no_of_divisions];
	}
}

// same as the American call

double american_put_option(int k, int i, double current_stock_price)
{
	double temp;
	if (momi2[k][i+no_of_divisions]!=-1) return momi2[k][i+no_of_divisions];
	if (k == no_of_divisions)
		return max(0.0, (strike_price - current_stock_price));
	else
	{
		if(momi2[k+1][i+no_of_divisions]==-1)
		{
			momi2[k+1][i+no_of_divisions]=american_put_option(k+1, i,current_stock_price);
		}
		if(momi2[k+1][i+1+no_of_divisions]==-1)
		{
			momi2[k+1][i+1+no_of_divisions]=american_put_option(k+1, i+1,current_stock_price*up_factor);
		}
		if(momi2[k+1][i-1+no_of_divisions]==-1)
		{
			momi2[k+1][i-1+no_of_divisions]=american_put_option(k+1, i-1,current_stock_price/up_factor);
		}
		temp=(uptick_prob*momi2[k+1][i+1+no_of_divisions]+downtick_prob*momi2[k+1][i-1+no_of_divisions]+ notick_prob*momi2[k+1][i+no_of_divisions])/R;
		momi2[k][i+no_of_divisions]=max(temp, (-current_stock_price + strike_price));
		return momi2[k][i+no_of_divisions];
	
	}
}

// written by professor R.S

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a=fabs(z);
	double t = 1.0/(1.0+a*p);
	double b = c2*exp((-z)*(z/2.0));
	double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
	n = 1.0-b*n;
	if ( z < 0.0 ) n = 1.0 - n;
	return n;
};

double option_price_put_black_scholes(const double& S,      // spot price
									  const double& K,      // Strike (exercise) price,
									  const double& r,      // interest rate
									  const double& sigma,  // volatility
									  const double& time){
	double time_sqrt = sqrt(time);
	double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1-(sigma*time_sqrt);
	return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
									   const double& K,       // strike (exercise) price,
									   const double& r,       // interest rate
									   const double& sigma,   // volatility
									   const double& time) {  // time to maturity
	double time_sqrt = sqrt(time);
	double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
	double d2 = d1-(sigma*time_sqrt);
	return S*N(d1) - K*exp(-r*time)*N(d2);
};

// same like the American call and put, just need to get rid of the max function
double european_call_option(int k, int i)
{
	if (momi1[k][i+no_of_divisions]!=-1) return momi1[k][i+no_of_divisions];
	if (k == no_of_divisions)
		return max(0.0, (initial_stock_price*pow(up_factor, ((float) i))) - strike_price);
	else
	{
		if(momi1[k+1][i+no_of_divisions]==-1)
		{
			momi1[k+1][i+no_of_divisions]=european_call_option(k+1, i);
		}
		if(momi1[k+1][i+1+no_of_divisions]==-1)
		{
			momi1[k+1][i+1+no_of_divisions]=european_call_option(k+1, i+1);
		}
		if(momi1[k+1][i-1+no_of_divisions]==-1)
		{
			momi1[k+1][i-1+no_of_divisions]=european_call_option(k+1, i-1);
		}
		momi1[k][i+no_of_divisions]=(uptick_prob*momi1[k+1][i+1+no_of_divisions]+downtick_prob*momi1[k+1][i-1+no_of_divisions]+ notick_prob*momi1[k+1][i+no_of_divisions])/R;
		return momi1[k][i+no_of_divisions];
	}
}

double european_put_option(int k, int i)
{
	if (momi2[k][i+no_of_divisions]!=-1) return momi2[k][i+no_of_divisions];
	if (k == no_of_divisions)
		return max(0.0, strike_price - (initial_stock_price*pow(up_factor, ((float) i))));
	else
	{
		if(momi2[k+1][i+no_of_divisions]==-1)
		{
			momi2[k+1][i+no_of_divisions]=european_put_option(k+1, i);
		}
		if(momi2[k+1][i+1+no_of_divisions]==-1)
		{
			momi2[k+1][i+1+no_of_divisions]=european_put_option(k+1, i+1);
		}
		if(momi2[k+1][i-1+no_of_divisions]==-1)
		{
			momi2[k+1][i-1+no_of_divisions]=european_put_option(k+1, i-1);
		}
		momi2[k][i+no_of_divisions]=(uptick_prob*momi2[k+1][i+1+no_of_divisions]+downtick_prob*momi2[k+1][i-1+no_of_divisions]+ notick_prob*momi2[k+1][i+no_of_divisions])/R;
		return momi2[k][i+no_of_divisions];
	}
}


// written by professor R.S
int main (int argc, char* argv[])
{
	
	int count;
	
	sscanf (argv[1], "%lf", &expiration_time);
	sscanf (argv[2], "%d", &no_of_divisions);
	sscanf (argv[3], "%lf", &risk_free_rate);
	sscanf (argv[4], "%lf", &volatility);
	sscanf (argv[5], "%lf", &initial_stock_price);
	sscanf (argv[6], "%lf", &strike_price);
	count=2*no_of_divisions;

// initialize my array
	for(int i=0;i<count;i++)
		for(int j=0;j<count;j++)
		{
			momi1[i][j]=-1;
			momi2[i][j]=-1;
		}
	
	R = exp(risk_free_rate*expiration_time/((float) no_of_divisions));
	up_factor = exp(volatility*sqrt((2*expiration_time)/((float) no_of_divisions)));
	uptick_prob = pow((sqrt(R) - 1/sqrt(up_factor))/(sqrt(up_factor)-1/sqrt(up_factor)),2);
	downtick_prob = pow((sqrt(up_factor) - sqrt(R))/(sqrt(up_factor)-1/sqrt(up_factor)),2);
	notick_prob = 1 - uptick_prob - downtick_prob;
	
	cout << "Recursive Trinomial American Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "R = " << R << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "Downtick Probability = " << downtick_prob << endl;
	cout << "Notick Probability = " << notick_prob << endl;
	cout << "--------------------------------------" << endl;
	cout<<endl;
	cout << "--------------------------------------" << endl;
	cout<<"First we compute American call and put"<<endl;
	double call_price1=american_call_option(0,0,initial_stock_price);
	double put_price1=american_put_option(0,0,initial_stock_price);
	cout << "Trinomial Price of an American Call Option = " << call_price1 << endl;
	cout << "Trinomial Price of an American Put Option = " << put_price1 << endl;
	cout << "--------------------------------------" << endl;
	cout << "Verifying Put-Call Parity: S+P-C = Kexp(-r*T)" << endl;
	cout <<  initial_stock_price << " + " << put_price1 << " - " << call_price1;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price1 - call_price1 << " ?=? " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	cout<<endl;cout<<endl;
	
	// initialize my array again
	for(int i=0;i<count;i++)
		for(int j=0;j<count;j++)
		{
			momi1[i][j]=-1;
			momi2[i][j]=-1;
		}
	
	cout<<"Then we compute European call and put"<<endl;
	
	// for european call and put, written by professor R.S
	cout << "--------------------------------------" << endl;
	double call_price = european_call_option(0, 0);
	cout << "Binomial Price of an European Call Option = " << call_price << endl;
	cout << "Call Price according to Black-Scholes = " <<
	option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
									volatility, expiration_time) << endl;
	double put_price = european_put_option(0, 0);
	cout << "Binomial Price of an European Put Option = " << put_price << endl;
	cout << "Put Price according to Black-Scholes = " <<
	option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
								   volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	cout << "Verifying Put-Call Parity: S+P-C = Kexp(-r*T)" << endl;
	cout <<  initial_stock_price << " + " << put_price << " - " << call_price;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price - call_price << " = " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	
}
