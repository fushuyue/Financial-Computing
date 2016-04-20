//
//  generateRV.hpp
//  IE525Project
//
//  Created by fushuyue on 4/19/16.
//  Copyright © 2016 Monte Carlo Simulation. All rights reserved.
//

#ifndef generateRV_h
#define generateRV_h

double min(double a,double b);

double max(double a,double b);

//double squaresum=0;
//double totalsum=0;

// // simulate uniform random variables

int Random();

double Uniform();

double gaussian_box_muller();

double cdf(double x);

double inverse(double u);

double Normalcdf(const double& z) ;

double callOptionPrice(double S,double t,double X,double r,double sigma,double q);


double control_monte_carlo_call_price(const int& n, const double& S, const double& K, const double& r, const double& q, const double& v, const double& T, double m);


double **sobol_points(unsigned N, unsigned D, char *dir_file);

double Asian_call_qmc(const int& n, const double& S, const double& K, const double& r, const double& q, const double& v, const double& T, double m,double l);



#endif /* generateRV_hpp */
