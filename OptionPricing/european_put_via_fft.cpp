// Pricing an European Call Option using Fast Fourier Transform (FFT)
// The technical details can be found in P. Carr and D. Madan, 
// "Option Valuation Using the Fast Fourier Transform," Journal of 
// Computational Finance, Volume 2, Number 4, Summer, 1999.
// This paper is available on the Compass Website.  It require some
// mathematical heavy-lifting.  Come see me if you have questions,
// as the nitty-gritty of the paper are a little beyond the scope of
// the class. The method, however, is not...
// 
// The grid size must be of the form 2^k for some integer k for the 
// FFT algorithm to work efficiently.   You will need to link with the
// Newmat library for the FFT algorithm etc.
// 

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include "newmat.h"
#include "newmatap.h"  // Need this for the FFT algorithm
#include "normdist.h"  // Need this for Odegaard's Black-Scholes routine
#define PI 3.141592654 // Declaring the constant PI for future use

using namespace std;

float risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility;
int grid_size;

float max(float a, float b) {
	return (b < a )? a:b;
}
// Prof. Odegaard's implementation of the Black-Scholes Call Price Formula
double option_price_call_black_scholes(const double& S,       // spot (underlying) price
									   const double& K,       // strike (exercise) price,
									   const double& r,       // interest rate
									   const double& sigma,   // volatility 
									   const double& time)	  // time to maturity 
{  
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt; 
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
};

// Odegaard's Black-Scholes Put C++ Code
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

// My implemenation of the Carr-Madan approach
double call_price_from_fft()
{
	double log_of_initial_stock_price = log(initial_stock_price);
	double log_of_strike_price = log(strike_price);
	
	// equation 20 of Carr & Madan's paper
	double b = log_of_strike_price;
	double lambda = (2*b)/((double) grid_size);
	
	// equation 23 of Carr & Madan's paper
	double eta = (2*PI)/(lambda * ((double) grid_size));
	
	// we are going to use an alpha of 2
	double alpha = 5.0;
	
	// characteristic function 6 and the expression for phi_T in my notes
	// for this we have to create a grid and do the necessary incantations
	complex <double> final_series[grid_size];
	
	for (int i = 0; i < grid_size; i++) {
		
		complex <double> numerator;
		complex <double> denominator;
		complex <double> phi, temp1, temp2;
		
		phi.real() = ((double) i)*eta;
		phi.imag() = -(alpha + 1);
		
		// Simpson's Rule weights
		double weight = 3.0 + pow(-1.0, (double)i+1);
		if (i == 0)
			weight -= 1.0;
		weight = weight/3.0;

		temp1.real() = (volatility*volatility*expiration_time)/2;
		temp1.imag() = 0.0;
		temp1 = temp1*phi*phi;
		
		temp2.real() = 0.0;
		temp2.imag() = log_of_initial_stock_price + (risk_free_rate - (volatility*volatility*0.5))*expiration_time;
		temp2 = temp2*phi;
		
		temp1 = exp(temp2 - temp1);
		
		numerator = weight*exp(-risk_free_rate*expiration_time)*temp1;
		
		denominator.real() = (alpha*alpha) + alpha - (((double) i*i)*eta*eta);
		denominator.imag() = (2*alpha + 1)*(((double) i)*eta);
		
		temp1.real() = 0.0;
		temp1.imag() = log_of_strike_price * eta * ((double) i);
		temp1 = exp(temp1);
		
		final_series[i] = (eta*temp1*numerator)/denominator;
	}
	
	// have to get the real and imaginary parts for FFT
	ColumnVector real_part(grid_size), imaginary_part(grid_size);
	for (int i = 1; i <= grid_size; i++) {
		real_part(i) = final_series[i-1].real();
		imaginary_part(i) = final_series[i-1].imag();
	}
	ColumnVector fft_real_part(grid_size), fft_imaginary_part(grid_size);
	FFT(real_part, imaginary_part, fft_real_part, fft_imaginary_part);
	return(fft_real_part(1)*exp(-alpha*log_of_strike_price)/PI);
}

// My implemenation of the Carr-Madan approach
double put_price_from_fft()
{
	double log_of_initial_stock_price = log(initial_stock_price);
	double log_of_strike_price = log(strike_price);
	
	// equation 20 of Carr & Madan's paper
	double b = log_of_strike_price;
	double lambda = (2*b)/((double) grid_size);
	
	// equation 23 of Carr & Madan's paper
	double eta = (2*PI)/(lambda * ((double) grid_size));
	
	// we are going to use an alpha of 2
	double alpha = -5.0;
	
	// characteristic function 6 and the expression for phi_T in my notes
	// for this we have to create a grid and do the necessary incantations
	complex <double> final_series[grid_size];
	
	for (int i = 0; i < grid_size; i++) {
		
		complex <double> numerator;
		complex <double> denominator;
		complex <double> phi, temp1, temp2;
		
		phi.real() = ((double) i)*eta;
		phi.imag() = -(alpha + 1);
		
		// Simpson's Rule weights
		double weight = 3.0 + pow(-1.0, (double)i+1);
		if (i == 0)
			weight -= 1.0;
		weight = weight/3.0;
		
		temp1.real() = (volatility*volatility*expiration_time)/2;
		temp1.imag() = 0.0;
		temp1 = temp1*phi*phi;
		
		temp2.real() = 0.0;
		temp2.imag() = log_of_initial_stock_price + (risk_free_rate - (volatility*volatility*0.5))*expiration_time;
		temp2 = temp2*phi;
		
		temp1 = exp(temp2 - temp1);
		
		numerator = weight*exp(-risk_free_rate*expiration_time)*temp1;
		
		denominator.real() = (alpha*alpha) + alpha - (((double) i*i)*eta*eta);
		denominator.imag() = (2*alpha + 1)*(((double) i)*eta);
		
		temp1.real() = 0.0;
		temp1.imag() = log_of_strike_price * eta * ((double) i);
		temp1 = exp(temp1);
		
		final_series[i] = (eta*temp1*numerator)/denominator;
	}
	
	// have to get the real and imaginary parts for FFT
	ColumnVector real_part(grid_size), imaginary_part(grid_size);
	for (int i = 1; i <= grid_size; i++) {
		real_part(i) = final_series[i-1].real();
		imaginary_part(i) = final_series[i-1].imag();
	}
	ColumnVector fft_real_part(grid_size), fft_imaginary_part(grid_size);
	FFT(real_part, imaginary_part, fft_real_part, fft_imaginary_part);
	return(fft_real_part(1)*exp(-alpha*log_of_strike_price)/PI);
}

																	   
int main (int argc, char* argv[])
{
	sscanf (argv[1], "%f", &expiration_time);
	sscanf (argv[2], "%d", &grid_size);
	sscanf (argv[3], "%f", &risk_free_rate);
	sscanf (argv[4], "%f", &volatility);
	sscanf (argv[5], "%f", &initial_stock_price);
	sscanf (argv[6], "%f", &strike_price);
	
	// checking if the grid size is 2^k for some k
	int temp = grid_size;
	while((temp >= 2) && (temp%2 == 0))
		temp = temp/2;

	if (temp == 1) {
		cout << "European Put Option Pricing by FFT" << endl;
		cout << "Expiration Time (Years) = " << expiration_time << endl;
		cout << "Grid Size = " << grid_size << endl;
		cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
		cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
		cout << "Initial Stock Price = " << initial_stock_price << endl;
		cout << "Strike Price = " << strike_price << endl;
		cout << "--------------------------------------" << endl;
		cout << "Price of an European Call Option Via Carr & Madan's FFT method = " << 
		call_price_from_fft() << endl;
		cout << "Price of an European Call from the Black-Scholes Formula = " << 
		option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time) << endl;
		cout << "--------------------------------------" << endl;
		cout << "Price of an European Put Option Via Carr & Madan's FFT method = " << 
		put_price_from_fft() << endl;
		cout << "Price of an European Put from the Black-Scholes Formula = " << 
		option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time) << endl;
		cout << "--------------------------------------" << endl;
	}
	else {
		cout << "Your Grid Size is " << grid_size << endl;
		cout << "Retry with something that is 2^k for some k" << endl;
		cout << "Exiting..." << endl;
	}

}

