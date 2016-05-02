// Pricing the Asian American Call Option using the Hull-White Interpolation Method

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

double up_factor, uptick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R, h;
int no_of_divisions, highest_index;

double v[500][501][2000];

double max(double a, double b) {
	return (b < a )? a:b;
}

double american_asian_call_hull_white_interpolation()
{
	double x1, x2, y1, y2;
	double average_if_we_saw_an_uptick, average_if_we_saw_a_downtick;
	double v_if_uptick, v_if_downtick;
	
	// fill the final stage values
	for (int i = 0; i <= no_of_divisions; i++) 
		for (int j = 0; j < (2*highest_index+1); j++) 
			v[no_of_divisions][i][j] = max(0.0, initial_stock_price*pow(((double) 2.718281828459045),  h *((double) j-highest_index)) - strike_price);
	
	
	// filling the rest of the v[][][] matrix
	for (int k = no_of_divisions-1; k >= 0; k--) {
		
		for (int i = 0; i <= k; i++) {
			// for example, k = 6, i = 1, denotes (6,-6+(2*i)) = (6,-4) in the lattice
			// i.e. (k, 2*i - k) is equivalent state on the binomial lattice
			
			for (int j = 0; j < (2*highest_index+1); j++) {
				
				// compute the average if you saw an uptick; ie next state is (k+1,i+1)
				average_if_we_saw_an_uptick = (((double) k)*initial_stock_price*pow(((double) 2.718281828459045), h*((double) j-highest_index)) + 
											   initial_stock_price*pow(up_factor,2*i-k+1))/((double) k+1);
				
				for (int l = 0; l < (2*highest_index+1); l++) 
					if ( ((initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index))) <= average_if_we_saw_an_uptick) &&
						((initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index+1))) > average_if_we_saw_an_uptick) ) {
						x1 = initial_stock_price*pow(((double) 2.718281828459045), h*(l-highest_index));
						x2 = initial_stock_price*pow(((double) 2.718281828459045), h*(l-highest_index+1));
						// if there is an uptick then (k,i) changes to (k+1, i+1) which is [k+1][i+1][] in the array notation I have used 
						y1 = v[k+1][i+1][l];  
						y2 = v[k+1][i+1][l+1];
					}
				v_if_uptick = y1 + ((y2-y1)/(x2-x1))*(average_if_we_saw_an_uptick - x1);
				
				// compute the average if you saw a downtick; ie next state is (k+1,i-1)
				average_if_we_saw_a_downtick = (((double) k)*initial_stock_price*pow(((double) 2.718281828459045), h*((double) j-highest_index)) + 
												initial_stock_price*pow(up_factor,2*i-k-1))/((double) k+1);
				
				for (int l = 0; l < (2*highest_index+1); l++)
					if ( ((initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index))) <= average_if_we_saw_a_downtick) &&
						((initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index+1))) > average_if_we_saw_a_downtick) ) {
						x1 = initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index));
						x2 = initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index+1));
						// if there is a downtick then (k,i) changes to (k+1,i-1) which is [k+1][i][] in the array notation I have used
						y1 = v[k+1][i][l];
						y2 = v[k+1][i][l+1];
					}
				v_if_downtick = y1 + ((y2-y1)/(x2-x1))*(average_if_we_saw_a_downtick - x1);	
				
				v[k][i][j] = max(initial_stock_price*pow(((double) 2.718281828459045), h*((double) j-highest_index)) - strike_price,
								 (uptick_prob*v_if_uptick + (1-uptick_prob)*v_if_downtick)/R);
			}
		}
	}
	// return the value v[0][0][middle-entry] which is v[0][0][highest_index] in the array notation I used
	return v[0][0][highest_index];
}

double american_asian_put_hull_white_interpolation()
{
	double x1, x2, y1, y2;
	double average_if_we_saw_an_uptick, average_if_we_saw_a_downtick;
	double v_if_uptick, v_if_downtick;
	
	// fill the final stage values
	for (int i = 0; i <= no_of_divisions; i++) 
		for (int j = 0; j < (2*highest_index+1); j++) 
			v[no_of_divisions][i][j] = max(0.0, strike_price - initial_stock_price*pow(((double) 2.718281828459045),  h *((double) j-highest_index)));
	
	
	// filling the rest of the v[][][] matrix
	for (int k = no_of_divisions-1; k >= 0; k--) {
		
		for (int i = 0; i <= k; i++) {
			// for example, k = 6, i = 1, denotes (6,-6+(2*i)) = (6,-4) in the lattice
			// i.e. (k, 2*i - k) is equivalent state on the binomial lattice
			
			for (int j = 0; j < (2*highest_index+1); j++) {
				
				// compute the average if you saw an uptick; ie next state is (k+1,i+1)
				average_if_we_saw_an_uptick = (((double) k)*initial_stock_price*pow(((double) 2.718281828459045), h*((double) j-highest_index)) + 
											   initial_stock_price*pow(up_factor,2*i-k+1))/((double) k+1);
				
				for (int l = 0; l < (2*highest_index+1); l++) 
					if ( ((initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index))) <= average_if_we_saw_an_uptick) &&
						((initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index+1))) > average_if_we_saw_an_uptick) ) {
						x1 = initial_stock_price*pow(((double) 2.718281828459045), h*(l-highest_index));
						x2 = initial_stock_price*pow(((double) 2.718281828459045), h*(l-highest_index+1));
						// if there is an uptick then (k,i) changes to (k+1, i+1) which is [k+1][i+1][] in the array notation I have used 
						y1 = v[k+1][i+1][l];  
						y2 = v[k+1][i+1][l+1];
					}
				v_if_uptick = y1 + ((y2-y1)/(x2-x1))*(average_if_we_saw_an_uptick - x1);
				
				// compute the average if you saw a downtick; ie next state is (k+1,i-1)
				average_if_we_saw_a_downtick = (((double) k)*initial_stock_price*pow(((double) 2.718281828459045), h*((double) j-highest_index)) + 
												initial_stock_price*pow(up_factor,2*i-k-1))/((double) k+1);
				
				for (int l = 0; l < (2*highest_index+1); l++)
					if ( ((initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index))) <= average_if_we_saw_a_downtick) &&
						((initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index+1))) > average_if_we_saw_a_downtick) ) {
						x1 = initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index));
						x2 = initial_stock_price*pow(((double) 2.718281828459045), h*((double) l-highest_index+1));
						// if there is a downtick then (k,i) changes to (k+1,i-1) which is [k+1][i][] in the array notation I have used
						y1 = v[k+1][i][l];
						y2 = v[k+1][i][l+1];
					}
				v_if_downtick = y1 + ((y2-y1)/(x2-x1))*(average_if_we_saw_a_downtick - x1);	
				
				v[k][i][j] = max(strike_price - initial_stock_price*pow(((double) 2.718281828459045), h*((double) j-highest_index)),
								 (uptick_prob*v_if_uptick + (1-uptick_prob)*v_if_downtick)/R);
			}
		}
	}
	// return the value v[0][0][middle-entry] which is v[0][0][highest_index] in the array notation I used
	return v[0][0][highest_index];
}

int main (int argc, char* argv[])
{
	
	sscanf (argv[1], "%lf", &expiration_time);
	sscanf (argv[2], "%d", &no_of_divisions);
	sscanf (argv[3], "%lf", &risk_free_rate);
	sscanf (argv[4], "%lf", &volatility);
	sscanf (argv[5], "%lf", &initial_stock_price);
	sscanf (argv[6], "%lf", &strike_price);
	sscanf (argv[7], "%lf", &h);
	
	up_factor = exp(volatility*sqrt(expiration_time/((double) no_of_divisions)));
	R = exp(risk_free_rate*expiration_time/((double) no_of_divisions));
	uptick_prob = (R - (1/up_factor))/(up_factor-(1/up_factor));
	
	// round number
	highest_index = max(floor(-volatility*sqrt(expiration_time/((double) no_of_divisions))*(no_of_divisions/h)),
						ceil(volatility*sqrt(expiration_time/((double) no_of_divisions))*(no_of_divisions/h))) + 1;
	
	
	cout << "American-Asian Option Pricing using Hull-White Interpolation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "h = " << h << endl;
	cout << "--------------------------------------" << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "Highest Index = " << highest_index << endl;
	cout << "R = " << R << endl;
	cout << "--------------------------------------" << endl;
	cout << "American-Asian Call Option Price by Hull-White Interpolation = " << american_asian_call_hull_white_interpolation() << endl;
	cout << "American-Asian Put Option Price by Hull-White Interpolation = " << american_asian_put_hull_white_interpolation() << endl;
	
}

